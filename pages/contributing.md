# Contributing

This is a sphinx project, which generates documentation for Carleton College BIOL.338 Genomics & Bioinformatics.

**To contribute, you must:**

1. Be able to edit markdown (easy!)
2. Be able to run sphinx, to generate webpages from markdown (easy!)
3. Be able to use git to upload file changes to github (also easy!)

**How do these tools fit together?**

- We write the docs in markdown, then use sphinx to render beatutiful webpages from them.
- We use Git/Github to store the current version of the project online.
- We can also use [ReadTheDocs](https://readthedocs.org/) to automatically build/serve our website, based on the current files stored on Github!

### Step 1: Markdown

[Markdown](http://commonmark.org) is a fairly simple specification for formatting text.

Markdown lets you indicate that things should be bolded or italicized or seen as headers or links in "plain-text" (ie, by typing certain characters rather than pushing application-specific buttons).

You can even write $\LaTeX$ in Markdown, by enclosing your Latex stuff in `$`.

> For syntax reference: see [commonmark.org/help](http://commonmark.org/help/). (Though keep in mind some Markdown rendering "engines", like the one on GitHub, support slightly different syntax, beyond and even different than the syntax specified by CommonMark.)


#### Using Markdown

A nice way to write Markdown and see it rendered into a formatted documnt is to use the [Atom](https://atom.io/) text editor, with the fantastic extension [Markdown Preview Enhanced](https://shd101wyy.github.io/markdown-preview-enhanced/#/).



### Step 2: Sphinx

#### Installing Sphinx

1. Install sphinx: `pip install sphinx`
2. Install the markdown extension: `pip install recommonmark`
3. Install the ReadTheDocs Sphinx theme: `pip install sphinx_rtd_theme`

> Note: `pip` commands are run at the command line / in Terminal.

#### Using Sphinx

Prepare materials:
- Write your markdown documents and put them somewhere, like the `pages/` folder.
- To control which files are included: edit the file `index.rst`
- To alter other settings: edit the file `conf.py`

Build the website:
- Run `make clean` in terminal to empty out the `_build/` directory
- Run `make html` in terminal to generate html from the markdown files.
- Type `open _build/html/index.html` in terminal to view your website with your browser (or open that file using Finder.)

> Note: A good example project, for inspiration on using sphinx/markdown is the [Requests](https://github.com/requests/requests) package, which is on Github. Look in their `docs/` directory, find files like `index.rst` and press the `raw` button.


### Step 3: Github

When you're ready to push your changes to github.

#### With Permissions

If you have permissions to push to the Github repository, commiting and pushing your changes just requires a few lines:

```bash
# see what has changed
git status

# stage modified files
git add -u

# add any new files
git add <filename> <filename>

# commit changes
git commit -m "< write a msg like 'modify protocol 5'>"

# push changes to remote server (github)
git push
```

#### Without Permissions

If you don't have Github pushing permission for this project, you should fork the project, commit your changes, then open a pull request. This is the way open source software gets built. Hooray for open source!!!

Here's a good resource on this workflow:
https://gist.github.com/Chaser324/ce0505fbed06b947d962