[CmdletBinding()]
param([string]$ContainerName = 'colab-gpu')

$ErrorActionPreference = 'Stop'
$dockerExe = 'C:\Program Files\Docker\Docker\resources\bin\docker.exe'
if (-not (Test-Path -LiteralPath $dockerExe)) { throw 'Docker Desktop is not installed.' }
if (& $dockerExe ps --format '{{.Names}}' | Where-Object { $_ -eq $ContainerName }) {
    & $dockerExe stop $ContainerName | Out-Null
    Write-Host "Stopped $ContainerName. C:\colab-work is preserved."
} else { Write-Host "$ContainerName is already stopped." }
