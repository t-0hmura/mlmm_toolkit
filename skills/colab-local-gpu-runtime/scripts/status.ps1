[CmdletBinding()]
param([string]$ContainerName = 'colab-gpu', [int]$HostPort = 9000)

$ErrorActionPreference = 'Stop'
$dockerExe = 'C:\Program Files\Docker\Docker\resources\bin\docker.exe'
if (-not (Test-Path -LiteralPath $dockerExe)) { throw 'Docker Desktop is not installed.' }
Write-Host 'Host GPU:'
& nvidia-smi.exe --query-gpu=name,driver_version,memory.total --format=csv,noheader
Write-Host 'Container:'
& $dockerExe ps --filter "name=$ContainerName" --format 'table {{.Names}}\t{{.Status}}\t{{.Ports}}'
if (-not (& $dockerExe ps --format '{{.Names}}' | Where-Object { $_ -eq $ContainerName })) { throw "Container '$ContainerName' is not running. Run start.ps1." }
Write-Host 'Container GPU:'
& $dockerExe exec $ContainerName nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
$servers = (& $dockerExe exec $ContainerName bash -lc 'jupyter server list') -join "`n"
$match = [regex]::Matches($servers, 'http://127\.0\.0\.1:9000/\?token=[A-Za-z0-9]+') | Select-Object -Last 1
if (-not $match) { throw 'Jupyter token URL is not ready. Wait and rerun status.ps1.' }
$url = $match.Value -replace '^http://127\.0\.0\.1:9000', "http://localhost:$HostPort"
Write-Host "`nColab local runtime URL:"
Write-Host $url -ForegroundColor Cyan
Write-Host "`nIn Colab: Connect -> Connect to a local runtime -> paste the URL above."
