$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$outputPath = [System.IO.Path]::GetFullPath((Join-Path $repoRoot '..\siwenyujuan-site'))
Push-Location $repoRoot
try {
    conda run --no-capture-output -n mkdocs python site-config/hooks.py build --strict --config-file site-config/mkdocs.yml
    if ($LASTEXITCODE -ne 0) {
        throw "MkDocs build failed with exit code $LASTEXITCODE."
    }
    Write-Host "Site generated at: $outputPath" -ForegroundColor Green
}
finally {
    Pop-Location
}
