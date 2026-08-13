# Colab connection and troubleshooting

## Connect

1. Run `powershell -ExecutionPolicy Bypass -File scripts/start.ps1` and copy the token URL.
2. Open a trusted pdb2reaction or mlmm-toolkit notebook in Colab.
3. Select **Connect ▸ Connect to a local runtime**.
4. Paste the URL, connect, run **Installation**, then run **Launch GUI**.

The UI remains in Colab. Python, files, and GPU computation run locally. `/content` maps to `C:\colab-work`.

## Installation modes

- Published version: retain the version tag and install the pinned release.
- Debug bundle: set the version field to `debug` and keep the matching source ZIP in `C:\colab-work`.
- Use ORB for the simplest no-login check.

## Troubleshoot

- **Connect stays disabled:** include the complete `?token=...` URL. Use Chrome or Edge because some embedded browsers block `localhost`.
- **No GPU:** run `status.ps1`; host and container must name the same NVIDIA GPU.
- **Docker will not start after first installation:** restart Windows once, open Docker Desktop, and run `start.ps1`.
- **Release not found on PyPI:** use the matching debug ZIP and set the version field to `debug`.
- **Files disappear:** save under `C:\colab-work`, not elsewhere in the container.
- **Drive mount fails:** local runtimes do not support `drive.mount()`.

Never publish the Jupyter token. Execute only trusted notebooks.
