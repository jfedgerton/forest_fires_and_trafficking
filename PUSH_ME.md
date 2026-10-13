# How to push this to GitHub

This archive is a complete git repository with its history. Extract it, point it
at the remote, push.

```bash
tar xzf forest_fires_and_trafficking.tar.gz
cd forest_fires_and_trafficking
git remote add origin https://github.com/jfedgerton/forest_fires_and_trafficking.git
git branch -M main
git push -u origin main
```

If the remote already has a commit (a README created with the repo), either

```bash
git pull --rebase origin main && git push -u origin main
```

or, since nothing here depends on that README,

```bash
git push -u origin main --force
```

Three commits carry the full recovery history. `MANIFEST.md` is the thing to
read first — it marks which files are byte-identical from Dropbox, which are
reconstructed and have never run, and which are gone.

No data is included. `.gitignore` excludes it deliberately.

You can delete this file after pushing.
