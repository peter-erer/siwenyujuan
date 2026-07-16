$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $repoRoot
try {
    Write-Host 'Starting site at http://127.0.0.1:8000/siwenyujuan/' -ForegroundColor Cyan
    conda run --no-capture-output -n mkdocs python site-config/hooks.py serve --config-file site-config/mkdocs.yml
}
finally {
    Pop-Location
}
