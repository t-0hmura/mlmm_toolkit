[CmdletBinding()]
param(
    [string]$ContainerName = 'colab-gpu',
    [int]$HostPort = 9000,
    [string]$WorkDir = 'C:\colab-work',
    [string]$BundleRoot = ''
)

$ErrorActionPreference = 'Stop'
if (-not $BundleRoot) { $BundleRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path }
$dockerExe = 'C:\Program Files\Docker\Docker\resources\bin\docker.exe'
$dockerDesktop = 'C:\Program Files\Docker\Docker\Docker Desktop.exe'
if (-not (Test-Path -LiteralPath $dockerExe)) { throw 'Docker Desktop is not installed. Run setup.ps1 first.' }
if (-not (Get-Process 'Docker Desktop' -ErrorAction SilentlyContinue)) { Start-Process -FilePath $dockerDesktop -WindowStyle Hidden }
$deadline = (Get-Date).AddMinutes(3)
do {
    $previous = $ErrorActionPreference; $ErrorActionPreference = 'SilentlyContinue'
    & $dockerExe info --format '{{.ServerVersion}}' *> $null; $dockerReady = $LASTEXITCODE
    $ErrorActionPreference = $previous
    if ($dockerReady -eq 0) { break }
    Start-Sleep -Seconds 5
} while ((Get-Date) -lt $deadline)
if ($dockerReady -ne 0) { throw 'Docker did not become ready.' }

New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
foreach ($name in @('pdb2reaction_colab.ipynb','mlmm_colab.ipynb','pdb2reaction-src.zip','mlmm-toolkit-src.zip','pdb2reaction_inputs','mlmm_inputs')) {
    $source = Join-Path $BundleRoot $name
    if (Test-Path -LiteralPath $source) { Copy-Item -LiteralPath $source -Destination $WorkDir -Recurse -Force }
}
if (-not (& $dockerExe ps -a --format '{{.Names}}' | Where-Object { $_ -eq $ContainerName })) { throw "Container '$ContainerName' is missing. Run setup.ps1 first." }
& $dockerExe start $ContainerName | Out-Null
# A restarted container can briefly report the token from the previous Jupyter
# process before its entrypoint finishes replacing that server.  Wait for the
# reported URL to remain unchanged across several probes so the URL printed
# below is the one Colab can actually use.
Start-Sleep -Seconds 10
$tokenDeadline = (Get-Date).AddMinutes(2)
$stableToken = $null
$stableCount = 0
do {
    $previous = $ErrorActionPreference; $ErrorActionPreference = 'SilentlyContinue'
    $servers = (& $dockerExe exec $ContainerName bash -lc 'jupyter server list' 2>$null) -join "`n"
    $probe = [regex]::Matches($servers, 'http://127\.0\.0\.1:9000/\?token=[A-Za-z0-9]+') | Select-Object -Last 1
    $ErrorActionPreference = $previous
    if ($probe -and $probe.Value -eq $stableToken) {
        $stableCount += 1
    } elseif ($probe) {
        $stableToken = $probe.Value
        $stableCount = 1
    } else {
        $stableToken = $null
        $stableCount = 0
    }
    if ($stableCount -ge 3) { break }
    Start-Sleep -Seconds 2
} while ((Get-Date) -lt $tokenDeadline)
if ($stableCount -lt 3) { throw 'Jupyter token URL did not become stable.' }
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'status.ps1') -ContainerName $ContainerName -HostPort $HostPort
