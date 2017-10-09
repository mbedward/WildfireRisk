# Private (non-exported) functions

# Takes a vector x and returns a vector of length n
# by recycling or truncating as required.
#
.fixed_length_vector <- function(x, n) {
  nx <- length(x)

  if (nx == 0) {
    xname <- deparse(substitute(x))
    stop("Argument ", xname, " is required")
  }

  if (nx < n)
    rep(x, length.out = n)
  else
    x[1:n]
}


# Checks that a function has a given number of arguments.
#
.check_nargs <- function(f, n) {
  nf <- length(formals(f))
  if (nf != n) {
    fname <- deparse(substitute(f))
    stop("Function ", fname, " takes ", nf, " arguments but should take ", n)
  }
}


# Retrieve an 'sf' object of vector data.
# If obj is an 'sf' object, just return it;
# if it is a string, treat it as the path to a file to read.
#
.get_sf_object <- function(obj) {
  if (is.character(obj))
    obj <- st_read(obj[1], quiet = TRUE)

  if (!inherits(obj, "sf")) {
    objname <- deparse(substitute(obj))
    stop("Argument ", objname, " must be either an sf object\n",
         "  or the path to a file of vector data to read")
  }

  obj
}

