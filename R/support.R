# support.R -- a miscellaneous collection of supporting functions

# Copyright (C) 2016 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

ctext <- function (str, n) {
    # Centers the string str in a field of n blanks
    if (nchar(str) == n) return(str)
    if (nchar(str) > n) return(substr(str, 1, n))
    
    nleft <- floor((n - nchar(str))/2)
    nright <- n - (nchar(str) + nleft)
    left_pad <- paste(rep(" ", nleft), collapse="")
    right_pad <- paste(rep(" ", nright), collapse="")
    pstr <- paste(left_pad, str, right_pad, sep="")
    pstr
}

printf <- function (...) { cat(sprintf(...)); }
println <- function (...) { cat(sprintf(...)); cat("\n"); }

