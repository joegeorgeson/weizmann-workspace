### Here are some tricks for using python code from R

How to setup environment - note that users can specify the python (and therefore environment)
```
library(reticulate)

path_to_python <- "~/.virtualenvs/r-reticulate/bin/python"
use_python(path_to_python)
```

Now you can use `pandas` to import a `pickle`
```
pd <- import("pandas")

pickle.file <- "~/my_pickle_file.pd"
x <- pd$read_pickle(pickle.file)
```
