# Contributing

This is a sphinx project, which generates documentation for Carleton College BIOL.338 Genomics & Bioinformatics.

To contribute, you must:

1. Be able to edit markdown (easy!)
2. Be able to run sphinx, to generate webpages from markdown (easy!)
3. Be able to use git to upload file changes to github (also easy!)

### Step 1: Markdown

**Learning:**
- [Markdown](http://commonmark.org) is a simple way to format text.

- You can indicate in plain-text (ie, by typing special symbols rather than pushing application-specific buttons) that things should be bolded or italicized or seen as headers or links.

- See [commonmark.org/help](http://commonmark.org/help/) for quick reference. (Though note some Markdown "engines" support other syntax, beyond and even different than the syntax specified by CommonMark.)

**Using:**

- A nice way to edit markdown and see it rendered into HTML is to use the [Atom](https://atom.io/) text editor, with the fantastic extension [Markdown Preview Enhanced](https://shd101wyy.github.io/markdown-preview-enhanced/#/).

### Step 2: Sphinx

**Installing:**
1. Install sphinx: `pip install sphinx`
2. Install the markdown extension: `pip install recommonmark`
3. Install the ReadTheDocs Sphinx theme: `pip install sphinx_rtd_theme`

**Using:**
- Edit the file `index.rst` to alter which files are included
- Edit the file `conf.py` to alter settings
- Run `make clean` in terminal to empty out the `_build/` directory
- Run `make html` in terminal to generate html from the markdown files.
- Type `open _build/html/index.html` in terminal to view your website with your browser.

> Note: A good example project, for inspiration on using sphinx/markdown is the [Requests](https://github.com/requests/requests) package, which is on Github. Look in their `docs/` directory, find files like `index.rst` and press the `raw` button.

### Step 3: github

When you're ready to push your changes to github:

```bash
# see what has changed
git status

# stage modified files
git add -u

# add any new files
git add <filename> <filename>

# commit changes
git commit -m "<msg like 'modify to protocol 5'>"

# push changes to remote server (github)
git push
```
