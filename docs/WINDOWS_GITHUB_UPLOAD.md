# Uploading PfVarX to GitHub from Windows 11

Target repository: `https://github.com/hadish-bekuretsion/PfVarX_Pipeline`
Local project path: `C:\Users\woldu\Downloads\Hadish\PfVarX_Pipeline`
Identity: `hadish-bekuretsion <hadish2028@gmail.com>`

---

## 0. Check that the folder is not cloud-synced

Your path is `C:\Users\woldu\Downloads\Hadish\PfVarX_Pipeline`, which is a local
folder — good. Git repositories and file-sync services (OneDrive, Dropbox,
Google Drive) do not coexist well, because a repository contains thousands of
small files under `.git\` that change constantly during a commit, and the
syncer can lock `.git\index` mid-write:

- `fatal: unable to write new index file`
- conflict copies such as `.git\index (1)` that git cannot read
- a repository corrupted by a commit interrupted mid-sync

Some Windows installations redirect the Desktop into OneDrive. Confirm yours
does not:

```powershell
echo $env:USERPROFILE
[Environment]::GetFolderPath("Desktop")
```

If the second command prints a path containing `OneDrive`, the folder is synced
and you should work elsewhere, for example:

```powershell
mkdir C:\dev
robocopy "C:\Users\woldu\Downloads\Hadish\PfVarX_Pipeline" "C:\dev\PfVarX_Pipeline" /E /XD .git __pycache__ .pytest_cache .ruff_cache
cd C:\dev\PfVarX_Pipeline
```

`robocopy` exits with code 1 when it copies files successfully; that is normal,
not an error. If the printed path is `C:\Users\woldu\Desktop`, you are fine and
can stay where you are.

---

## 1. Install Git for Windows

Check whether it is already there:

```powershell
git --version
```

If that fails, download and run the installer from
<https://git-scm.com/download/win>. Accept the defaults. The default editor is
Vim, which confuses many people; if you prefer, choose **Use Notepad as Git's
default editor** during setup. Everything else can stay default.

Close and reopen PowerShell afterwards so `git` is on the `PATH`.

---

## 2. Set your identity

These must match the account that will own the repository.

```powershell
git config --global user.name  "hadish-bekuretsion"
git config --global user.email "hadish2028@gmail.com"
git config --global init.defaultBranch main
```

Verify:

```powershell
git config --global --list
```

### Line endings

The repository ships a `.gitattributes` that stores text as LF. On Windows,
tell git not to second-guess it:

```powershell
git config --global core.autocrlf false
```

This is important. With the Windows default (`true`), every text file would be
rewritten to CRLF on checkout and git would report the entire tree as modified
even when you changed nothing.

---

## 3. Create the empty repository on GitHub

Open <https://github.com/new> in a browser.

- **Repository name:** `PfVarX_Pipeline`
- **Visibility:** Private is sensible while the software is unvalidated.
- **Do not** tick *Add a README file*, *Add .gitignore*, or *Choose a license*.

You already have all three files. If GitHub creates its own, the first push is
rejected with `Updates were rejected because the remote contains work that you
do not have locally`, and you then have to reconcile two unrelated histories.
An empty repository avoids that entirely.

Click **Create repository**. Leave the page open; you need the URL.

---

## 4. Create a Personal Access Token

GitHub no longer accepts account passwords for git operations. You need a
token.

1. Go to <https://github.com/settings/tokens>
2. **Generate new token → Generate new token (classic)**
3. **Note:** `PfVarX push from Windows`
4. **Expiration:** 90 days is reasonable
5. **Scopes:** tick **`repo`** (the whole group). Nothing else is required.
6. **Generate token**
7. Copy the token — it starts with `ghp_`. GitHub shows it once only.

Treat the token like a password. Do not paste it into a file inside the
project, and do not commit it.

---

## 5. Initialise the repository and commit

Run in PowerShell from the project root.

```powershell
cd C:\dev\PfVarX_Pipeline

git init -b main
git add -A
git status
```

**Read the `git status` output before committing.** You should see roughly 114
files. Confirm that none of the following appear:

- `__pycache__\` or `.pyc` files
- `.pytest_cache\` or `.ruff_cache\`
- `benchmarks\results\` (1.7 MB of generated output)
- any `.bam`, `.bai`, `.cram`, `.fastq.gz`, `.vcf.gz`, `.bcf`

The bundled `.gitignore` excludes all of these. If any still appear, stop and
check that `.gitignore` was copied across.

Then commit:

```powershell
git commit -m "PfVarX 0.9.0: audited SNP/small-INDEL pipeline for P. falciparum"
```

Confirm the author is correct:

```powershell
git log -1 --format="%an <%ae>"
```

Expected: `hadish-bekuretsion <hadish2028@gmail.com>`

---

## 6. Connect the remote and push

```powershell
git remote add origin https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
git push -u origin main
```

A window appears asking for credentials:

- **Username:** `hadish-bekuretsion`
- **Password:** paste the `ghp_...` token (right-click pastes in PowerShell;
  `Ctrl+V` may not work in the credential dialog)

If the push succeeds you will see something like:

```
 * [new branch]      main -> main
branch 'main' set up to track 'origin/main'.
```

### Avoid retyping the token every time

Windows ships with Git Credential Manager, which usually stores it
automatically. If it does not, enable it explicitly:

```powershell
git config --global credential.helper manager
```

---

## 7. Verify

```powershell
git remote -v
git log --oneline
git status
```

- `git remote -v` shows `origin` pointing at your repository for both fetch and push.
- `git log --oneline` shows your commit(s).
- `git status` reports `nothing to commit, working tree clean`.

Refresh <https://github.com/hadish-bekuretsion/PfVarX_Pipeline>. You should see
`README.md`, `src/`, `docs/`, `tests/`, `workflow/`, `config/`, `pyproject.toml`
and `LICENSE`.

---

## 8. If something goes wrong

**`remote origin already exists`**
```powershell
git remote set-url origin https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
```

**`Updates were rejected because the remote contains work that you do not have locally`**
You created the repository with a README or licence. Either delete the
repository on GitHub and recreate it empty, or merge the two histories:
```powershell
git pull origin main --allow-unrelated-histories
```
Resolve any conflicts, then `git push -u origin main`.

**`Authentication failed`**
The token is wrong, expired, or lacks the `repo` scope. Generate a new one.
To clear a cached bad credential, open *Credential Manager* in Windows, find
the entry under **Windows Credentials** beginning `git:https://github.com`,
and remove it.

**`fatal: unable to write new index file`**
OneDrive is holding a lock. Pause sync, or move the project out of OneDrive
(step 0).

**`LF will be replaced by CRLF` warnings**
`core.autocrlf` is still `true`. Run `git config --global core.autocrlf false`
and then `git add -A` again.

**Files appear as modified but you changed nothing**
Same cause as above. After setting `core.autocrlf false`, refresh the index:
```powershell
git rm --cached -r .
git reset --hard
```

---

## 9. Optional: preserve the original commit history

The project was already committed on a Linux machine with two commits. If you
want that history rather than a single fresh commit, use the bundle:

```powershell
git clone pfvarx.bundle PfVarX_Pipeline
cd PfVarX_Pipeline
git remote set-url origin https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
git push -u origin main
```

A bundle is a single file containing the full repository. It is safer than
copying a `.git` folder across filesystems.

---

## 10. Before you make it public

- `LICENSE` is present. Read it and confirm it is the licence you intend.
- `docs/PfVarX_Methodology.docx` and `.pdf` are tracked. They regenerate from
  `scripts/build_methodology_docx.sh`. If you would rather not carry build
  artefacts, add `docs/*.pdf` and `docs/*.docx` to `.gitignore` and run
  `git rm --cached docs/PfVarX_Methodology.docx docs/PfVarX_Methodology.pdf`
  before pushing.
- The audit report states plainly that all parameters are
  `PROVISIONAL_DEVELOPMENT` and the software is not clinically validated.
  Keep that statement visible; it is the honest description of the project's
  status.
