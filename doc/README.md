# README

The docs are available here as plain text files. Until they're made public,
you can view them in html by following the below steps:

1. Install [Sphinx](http://www.sphinx-doc.org/en/stable/tutorial.html).

```
$ pip install Sphinx
$ pip install sphinx_rtd_theme
```

2. Build the docs

```
$ git checkout docs
$ cd ${pipeline-directory}/doc
$ make html
```

3. Open the build in your web browser of choice

```
$ open _build/html/index.html
```
