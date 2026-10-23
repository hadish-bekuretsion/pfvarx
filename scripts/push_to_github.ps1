# push_to_github.ps1
#
# Prepare and push the PfVarX repository to GitHub from Windows.
#
# Usage (from PowerShell, in the project root):
#     .\scripts\push_to_github.ps1
#
# The script is deliberately cautious: it refuses to run if it would publish
# caches, generated output or genomic data, and it stops rather than guessing
# when something is ambiguous.

[CmdletBinding()]
param(
    [string]$RepoUrl = "https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git",
    [string]$UserName = "hadish-bekuretsion",
    [string]$UserEmail = "hadish2028@gmail.com",
    [string]$Branch = "main"
)

$ErrorActionPreference = "Stop"

function Fail($message) {
    Write-Host ""
    Write-Host "STOPPED: $message" -ForegroundColor Red
    exit 1
}

function Step($message) {
    Write-Host ""
    Write-Host "==> $message" -ForegroundColor Cyan
}

# ---------------------------------------------------------------------------
Step "Checking that git is installed"
# ---------------------------------------------------------------------------
try {
    $gitVersion = git --version
    Write-Host "    $gitVersion"
} catch {
    Fail "git is not installed or not on PATH. Install it from https://git-scm.com/download/win and reopen PowerShell."
}

# ---------------------------------------------------------------------------
Step "Checking the working directory"
# ---------------------------------------------------------------------------
if (-not (Test-Path "pyproject.toml") -or -not (Test-Path "src\pfvarx")) {
    Fail "This does not look like the PfVarX project root. Run the script from the folder that contains pyproject.toml and src\pfvarx."
}
$projectRoot = (Get-Location).Path
Write-Host "    $projectRoot"

$desktop = [Environment]::GetFolderPath("Desktop")
if ($desktop -match "OneDrive" -and $projectRoot -like "$desktop*") {
    Write-Host ""
    Write-Host "    WARNING: the project is inside a OneDrive-synced Desktop." -ForegroundColor Yellow
    Write-Host "    OneDrive can lock .git\index during a commit and corrupt the repository." -ForegroundColor Yellow
    $answer = Read-Host "    Continue anyway? (y/N)"
    if ($answer -ne "y") { Fail "Move the project to a local folder such as C:\dev\PfVarX_Pipeline and run the script again." }
}

# ---------------------------------------------------------------------------
Step "Configuring identity and line endings (repository-local)"
# ---------------------------------------------------------------------------
git config user.name  $UserName
git config user.email $UserEmail
# The repository ships .gitattributes that stores text as LF. On Windows the
# default core.autocrlf=true would rewrite every file to CRLF on checkout and
# report the whole tree as modified.
git config core.autocrlf false
Write-Host "    user.name  = $(git config user.name)"
Write-Host "    user.email = $(git config user.email)"
Write-Host "    core.autocrlf = $(git config core.autocrlf)"

# ---------------------------------------------------------------------------
Step "Initialising the repository"
# ---------------------------------------------------------------------------
if (Test-Path ".git") {
    Write-Host "    Repository already exists; reusing it."
} else {
    git init -b $Branch | Out-Null
    Write-Host "    Initialised on branch $Branch"
}

# ---------------------------------------------------------------------------
Step "Staging files"
# ---------------------------------------------------------------------------
git add -A
if ($LASTEXITCODE -ne 0) { Fail "git add failed." }

$staged = @(git diff --cached --name-only)
Write-Host "    $($staged.Count) files staged"

if ($staged.Count -eq 0) {
    Write-Host "    Nothing new to commit."
}

# ---------------------------------------------------------------------------
Step "Refusing to publish caches, output or genomic data"
# ---------------------------------------------------------------------------
$forbidden = @(
    '__pycache__',
    '\.pyc$',
    '\.pytest_cache',
    '\.ruff_cache',
    '\.egg-info',
    '^benchmarks/results/',
    '\.bam$',
    '\.bai$',
    '\.cram$',
    '\.crai$',
    '\.fastq$',
    '\.fastq\.gz$',
    '\.fq\.gz$',
    '\.vcf\.gz$',
    '\.bcf$'
)

$offenders = @()
foreach ($pattern in $forbidden) {
    $hits = @($staged | Where-Object { $_ -match $pattern })
    if ($hits.Count -gt 0) {
        $offenders += $hits
    }
}

if ($offenders.Count -gt 0) {
    Write-Host ""
    Write-Host "    These files must not be committed:" -ForegroundColor Red
    $offenders | Select-Object -Unique | ForEach-Object { Write-Host "      $_" -ForegroundColor Red }
    Fail "Unstage them with 'git reset', confirm .gitignore covers them, then run the script again."
}
Write-Host "    No caches, generated output or sequence data staged."

# ---------------------------------------------------------------------------
Step "Committing"
# ---------------------------------------------------------------------------
$pending = @(git diff --cached --name-only)
if ($pending.Count -gt 0) {
    $message = "PfVarX 0.9.0: audited SNP/small-INDEL pipeline for P. falciparum"
    git commit -m $message | Out-Null
    if ($LASTEXITCODE -ne 0) { Fail "git commit failed." }
    Write-Host "    Committed."
} else {
    Write-Host "    Nothing to commit."
}

Write-Host "    HEAD: $(git log -1 --format='%h %s')"
Write-Host "    Author: $(git log -1 --format='%an <%ae>')"

# ---------------------------------------------------------------------------
Step "Configuring the remote"
# ---------------------------------------------------------------------------
$existing = git remote 2>$null
if ($existing -contains "origin") {
    $current = git remote get-url origin
    if ($current -ne $RepoUrl) {
        Write-Host "    Updating origin from $current"
        git remote set-url origin $RepoUrl
    } else {
        Write-Host "    origin already points at $RepoUrl"
    }
} else {
    git remote add origin $RepoUrl
    Write-Host "    Added origin -> $RepoUrl"
}

# ---------------------------------------------------------------------------
Step "Pushing"
# ---------------------------------------------------------------------------
Write-Host "    You will be asked for credentials."
Write-Host "      Username: $UserName"
Write-Host "      Password: a Personal Access Token (ghp_...), NOT your account password."
Write-Host "      Create one at https://github.com/settings/tokens with the 'repo' scope."
Write-Host ""

git push -u origin $Branch
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "The push failed. Common causes:" -ForegroundColor Yellow
    Write-Host "  - The GitHub repository does not exist yet: create it empty at https://github.com/new" -ForegroundColor Yellow
    Write-Host "  - The repository was created with a README/licence, so the histories are unrelated." -ForegroundColor Yellow
    Write-Host "    If the remote holds only GitHub's auto-generated files, overwrite it:" -ForegroundColor Yellow
    Write-Host "        git push -u origin $Branch --force" -ForegroundColor Yellow
    Write-Host "    If the remote holds work you want to keep, merge instead:" -ForegroundColor Yellow
    Write-Host "        git pull origin $Branch --allow-unrelated-histories" -ForegroundColor Yellow
    Write-Host "        git push -u origin $Branch" -ForegroundColor Yellow
    Write-Host "  - The token is wrong, expired, or lacks the 'repo' scope." -ForegroundColor Yellow
    Write-Host "  - A cached bad credential: open Credential Manager, remove the git:https://github.com entry." -ForegroundColor Yellow
    exit 1
}

# ---------------------------------------------------------------------------
Step "Done"
# ---------------------------------------------------------------------------
Write-Host "    Pushed $Branch to $RepoUrl"
Write-Host "    Verify with:"
Write-Host "        git remote -v"
Write-Host "        git log --oneline"
Write-Host "        git status"
Write-Host ""
Write-Host "    Then open https://github.com/hadish-bekuretsion/PfVarX_Pipeline" -ForegroundColor Green
