# nullifier-pir

`nullifier-pir` serves PIR queries over fixed-width 32-byte nullifier snapshots.
The default snapshot shape is the 3317500 file:

```bash
cargo run -p nullifier-pir -- download \
  --url https://vote.fra1.cdn.digitaloceanspaces.com/snapshots/3317500/nullifiers.bin \
  --output data/nullifiers.bin
```

Start the local `ipir-sp` backend:

```bash
cargo run --release -p nullifier-pir -- serve \
  --snapshot-path data/nullifiers.bin \
  --backend local-ipir \
  --host 127.0.0.1 \
  --port 8080
```

The server exposes:

- `GET /health`
- `GET /meta`
- `POST /query` with backend-native query bytes

To compile the YPIR SimplePIR artifact backend pinned at commit `4f7ef3d`:

```bash
cargo check -p nullifier-pir --features ypir-artifact
```

## Packing Shape

The snapshot length is `1,597,627,296` bytes, or `49,925,853` nullifiers. The
crate packs `112` nullifiers into one SimplePIR item:

```text
112 nullifiers * 32 bytes * 8 bits = 28,672 bits
2048 coefficients * 14 bits       = 28,672 bits
```

The full snapshot therefore maps to `445,767` logical PIR items and pads to
`524,288` SimplePIR rows.

## Resource Notes

The local `ipir-sp` path stores the SimplePIR database in memory after encoding.
For the full snapshot this is roughly `524,288 * 2048 * sizeof(u16)`, about
2 GiB before preprocessing and Actix overhead. Preprocessing also allocates
InspiRING packing caches, so use a large-memory host for the full snapshot.
