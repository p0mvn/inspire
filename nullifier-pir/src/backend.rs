//! Backend abstraction for nullifier PIR servers.

use anyhow::{Context, Result};
use ipir_sp::client::IPIRClient;
use ipir_sp::params_for_simplepir;
use ipir_sp::server::{build_pack_preprocessed_blocks, IPIRServer};
use ipir_sp::YpirSchemeParams;
use serde::{Deserialize, Serialize};

use crate::encoding::ITEM_SIZE_BITS;
use crate::snapshot::NullifierSnapshot;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum BackendKind {
    LocalIpir,
    YpirArtifact,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendMetadata {
    pub backend: BackendKind,
    pub record_count: usize,
    pub pir_item_count: usize,
    pub db_rows: usize,
    pub db_cols: usize,
    pub item_size_bits: u64,
    pub setup_seed: u64,
}

pub trait PirBackend: Send + Sync {
    fn meta(&self) -> BackendMetadata;
    fn answer_query(&self, query: &[u8]) -> Result<Vec<u8>>;
}

pub fn seed_from_u64(value: u64) -> [u8; 32] {
    let mut seed = [0u8; 32];
    seed[..8].copy_from_slice(&value.to_le_bytes());
    seed
}

pub struct LocalIpirBackend {
    rlwe: &'static inspiring::RlweParams,
    ypir: YpirSchemeParams,
    record_count: usize,
    pir_item_count: usize,
    setup_seed: u64,
    server: IPIRServer<u16>,
    preprocessed: Vec<inspiring::PackPreprocessed<'static>>,
}

impl LocalIpirBackend {
    pub fn prepare(snapshot: &NullifierSnapshot, setup_seed: u64) -> Result<Self> {
        let (rlwe, ypir) = params_for_simplepir(snapshot.pir_row_count() as u64, ITEM_SIZE_BITS)
            .context("derive local ipir-sp SimplePIR parameters")?;
        Self::prepare_with_params(snapshot, setup_seed, rlwe, ypir)
    }

    pub fn prepare_with_params(
        snapshot: &NullifierSnapshot,
        setup_seed: u64,
        rlwe: inspiring::RlweParams,
        ypir: YpirSchemeParams,
    ) -> Result<Self> {
        let rlwe = Box::leak(Box::new(rlwe));
        let client = Box::leak(Box::new(IPIRClient::new(rlwe, &ypir)));
        let setup = client.generate_setup_simplepir_from_seed(seed_from_u64(setup_seed));
        let db = snapshot
            .coeff_iter(ypir.db_rows)
            .context("open snapshot coefficient iterator")?;
        let server = IPIRServer::<u16>::new(ypir.clone(), db, false, true);
        let offline = server.perform_offline_precomputation_simplepir(
            client.rlwe_params(),
            &setup.offline_query_polys,
        );
        let preprocessed = build_pack_preprocessed_blocks(
            client.rlwe_params(),
            &offline.crs_blocks,
            setup.key_pairs,
        )
        .context("build local ipir-sp preprocessing")?;

        Ok(Self {
            rlwe,
            ypir,
            record_count: snapshot.record_count(),
            pir_item_count: snapshot.pir_row_count(),
            setup_seed,
            server,
            preprocessed,
        })
    }
}

impl PirBackend for LocalIpirBackend {
    fn meta(&self) -> BackendMetadata {
        BackendMetadata {
            backend: BackendKind::LocalIpir,
            record_count: self.record_count,
            pir_item_count: self.pir_item_count,
            db_rows: self.ypir.db_rows,
            db_cols: self.ypir.db_cols,
            item_size_bits: self.ypir.item_size_bits,
            setup_seed: self.setup_seed,
        }
    }

    fn answer_query(&self, query: &[u8]) -> Result<Vec<u8>> {
        self.server
            .perform_full_online_computation_simplepir(self.rlwe, query, &self.preprocessed)
            .context("local ipir-sp query failed")
    }
}

pub enum Backend {
    Local(LocalIpirBackend),
    #[cfg(feature = "ypir-artifact")]
    YpirArtifact(ypir_artifact::YpirArtifactBackend),
}

impl Backend {
    pub fn prepare(
        kind: BackendKind,
        snapshot: &NullifierSnapshot,
        setup_seed: u64,
    ) -> Result<Self> {
        match kind {
            BackendKind::LocalIpir => Ok(Self::Local(LocalIpirBackend::prepare(
                snapshot, setup_seed,
            )?)),
            BackendKind::YpirArtifact => {
                #[cfg(feature = "ypir-artifact")]
                {
                    Ok(Self::YpirArtifact(
                        ypir_artifact::YpirArtifactBackend::prepare(snapshot, setup_seed)?,
                    ))
                }
                #[cfg(not(feature = "ypir-artifact"))]
                {
                    anyhow::bail!(
                        "backend ypir-artifact requires building with --features ypir-artifact"
                    )
                }
            }
        }
    }
}

impl PirBackend for Backend {
    fn meta(&self) -> BackendMetadata {
        match self {
            Self::Local(backend) => backend.meta(),
            #[cfg(feature = "ypir-artifact")]
            Self::YpirArtifact(backend) => backend.meta(),
        }
    }

    fn answer_query(&self, query: &[u8]) -> Result<Vec<u8>> {
        match self {
            Self::Local(backend) => backend.answer_query(query),
            #[cfg(feature = "ypir-artifact")]
            Self::YpirArtifact(backend) => backend.answer_query(query),
        }
    }
}

#[cfg(feature = "ypir-artifact")]
mod ypir_artifact {
    use super::*;
    use std::sync::Mutex;
    use ypir::params::{params_for_scenario_simplepir, DbRowsCols};
    use ypir::serialize::OfflinePrecomputedValues;
    use ypir::server::YServer;

    pub struct YpirArtifactBackend {
        params: &'static ypir_spiral::params::Params,
        record_count: usize,
        pir_item_count: usize,
        setup_seed: u64,
        server: YServer<'static, u16>,
        offline: Mutex<OfflinePrecomputedValues<'static>>,
    }

    impl YpirArtifactBackend {
        pub fn prepare(snapshot: &NullifierSnapshot, setup_seed: u64) -> Result<Self> {
            let params = Box::leak(Box::new(params_for_scenario_simplepir(
                snapshot.pir_row_count() as u64,
                ITEM_SIZE_BITS,
            )));
            let db = snapshot
                .coeff_iter(params.db_rows())
                .context("open snapshot coefficient iterator")?;
            let server = YServer::<u16>::new(params, db, false, true);
            let offline = server.perform_offline_precomputation_simplepir(None, None, None);

            Ok(Self {
                params,
                record_count: snapshot.record_count(),
                pir_item_count: snapshot.pir_row_count(),
                setup_seed,
                server,
                offline: Mutex::new(offline),
            })
        }
    }

    impl PirBackend for YpirArtifactBackend {
        fn meta(&self) -> BackendMetadata {
            BackendMetadata {
                backend: BackendKind::YpirArtifact,
                record_count: self.record_count,
                pir_item_count: self.pir_item_count,
                db_rows: self.params.db_rows(),
                db_cols: self.params.db_cols_simplepir(),
                item_size_bits: ITEM_SIZE_BITS,
                setup_seed: self.setup_seed,
            }
        }

        fn answer_query(&self, query: &[u8]) -> Result<Vec<u8>> {
            let offline = self
                .offline
                .lock()
                .map_err(|_| anyhow::anyhow!("YPIR offline cache mutex poisoned"))?;
            Ok(self
                .server
                .perform_full_online_computation_simplepir(&offline, query))
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use inspiring::{GadgetParams, RlweParams};
    use tempfile::NamedTempFile;

    use super::*;

    fn tiny_ypir(db_rows: usize, db_cols: usize) -> YpirSchemeParams {
        YpirSchemeParams {
            num_items: db_rows as u64,
            item_size_bits: (db_cols * 2) as u64,
            poly_len: 8,
            db_dim_1: 0,
            db_dim_2: 1,
            instances: db_cols / 8,
            db_rows,
            db_cols,
            p: 4,
            q_prime_1: 16,
            q_prime_2: 257,
            q2_bits: 8,
            t_exp_left: 3,
            t_exp_right: 2,
        }
    }

    #[test]
    fn local_backend_reports_snapshot_shape_with_tiny_params() {
        let mut file = NamedTempFile::new().expect("temp file");
        file.write_all(&[9u8; 32]).expect("write snapshot");
        let snapshot = NullifierSnapshot::open(file.path()).expect("open snapshot");
        let rlwe = RlweParams::new(
            8,
            12289,
            4,
            3.2,
            GadgetParams {
                bits_per: 3,
                ell: 5,
            },
        )
        .expect("valid params");

        let backend = LocalIpirBackend::prepare_with_params(&snapshot, 7, rlwe, tiny_ypir(8, 8))
            .expect("prepare backend");
        let meta = backend.meta();

        assert_eq!(meta.backend, BackendKind::LocalIpir);
        assert_eq!(meta.record_count, 1);
        assert_eq!(meta.pir_item_count, 1);
        assert_eq!(meta.db_rows, 8);
        assert_eq!(meta.db_cols, 8);
    }
}
