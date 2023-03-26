This template demonstrates many of the bells and whistles of the `reprex::reprex_document()` output format. The YAML sets many options to non-default values, such as using `#;-)` as the comment in front of output.

## Code style

Since `style` is `TRUE`, this difficult-to-read code (look at the `.Rmd` source file) will be restyled according to the Tidyverse style guide when it’s rendered. Whitespace rationing is not in effect!

``` r
x=1;y=2;z=x+y;z
#;-) [1] 3
```

## Quiet tidyverse

The tidyverse meta-package is quite chatty at startup, which can be very useful in exploratory, interactive work. It is often less useful in a reprex, so by default, we suppress this.

However, when `tidyverse_quiet` is `FALSE`, the rendered result will include a tidyverse startup message about package versions and function masking.

``` r
library(tidyverse)
#;-) ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
#;-) ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
#;-) ✔ tibble  3.1.8      ✔ dplyr   1.0.10
#;-) ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
#;-) ✔ readr   2.1.3      ✔ forcats 0.5.2 
#;-) ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#;-) ✖ dplyr::filter() masks stats::filter()
#;-) ✖ dplyr::lag()    masks stats::lag()
```

## Chunks in languages other than R

Remember: knitr supports many other languages than R, so you can reprex bits of code in Python, Ruby, Julia, C++, SQL, and more. Note that, in many cases, this still requires that you have the relevant external interpreter installed.

Let’s try Python!

``` python
x = 'hello, python world!'
print(x.split(' '))
#;-) ['hello,', 'python', 'world!']
```

And bash!

``` bash
echo "Hello Bash!";
pwd;
ls | head;
#;-) Hello Bash!
#;-) /home/georges/R/gell/R
#;-) drawings
#;-) inter_ells_reprex.Rmd
#;-) inter_ells_reprex_std_out_err.txt
#;-) inter_ells.Rmd
#;-) projective.R
```

Write a function in C++, use Rcpp to wrap it and …

``` cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
```

then immediately call your C++ function from R!

``` r
timesTwo(1:4)
#;-) [1] 2 4 6 8
```

## Standard output and error

Some output that you see in an interactive session is not actually captured by rmarkdown, when that same code is executed in the context of an `.Rmd` document. When `std_out_err` is `TRUE`, `reprex::reprex_render()` uses a feature of `callr:r()` to capture such output and then injects it into the rendered result.

Look for this output in a special section of the rendered document (and notice that it does not appear right here).

``` r
system2("echo", args = "Output that would normally be lost")
```

## Session info

Because `session_info` is `TRUE`, the rendered result includes session info, even though no such code is included here in the source document.

<details style="margin-bottom:10px;">
<summary>
Standard output and standard error
</summary>

``` sh
✖ Install the styler package in order to use `style = TRUE`.
running: python  -c "x = 'hello, python world!'
print(x.split(' '))"
running: bash  -c 'echo "Hello Bash!";
pwd;
ls | head;'
Building shared library for Rcpp code chunk...
Output that would normally be lost
```

</details>
<details style="margin-bottom:10px;">
<summary>
Session info
</summary>

``` r
sessioninfo::session_info()
#;-) ─ Session info ───────────────────────────────────────────────────────────────
#;-)  setting  value
#;-)  version  R version 4.2.2 Patched (2022-11-10 r83330)
#;-)  os       Ubuntu 22.04.1 LTS
#;-)  system   x86_64, linux-gnu
#;-)  ui       X11
#;-)  language en_CA:en
#;-)  collate  en_CA.UTF-8
#;-)  ctype    en_CA.UTF-8
#;-)  tz       America/New_York
#;-)  date     2023-01-20
#;-)  pandoc   2.19.2 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
#;-) 
#;-) ─ Packages ───────────────────────────────────────────────────────────────────
#;-)  package       * version date (UTC) lib source
#;-)  assertthat      0.2.1   2019-03-21 [1] CRAN (R 4.2.0)
#;-)  backports       1.4.1   2021-12-13 [1] CRAN (R 4.2.1)
#;-)  broom           1.0.1   2022-08-29 [1] CRAN (R 4.2.1)
#;-)  cellranger      1.1.0   2016-07-27 [1] CRAN (R 4.2.0)
#;-)  cli             3.4.1   2022-09-23 [1] CRAN (R 4.2.1)
#;-)  colorspace      2.0-3   2022-02-21 [1] CRAN (R 4.2.1)
#;-)  crayon          1.5.2   2022-09-29 [1] CRAN (R 4.2.1)
#;-)  DBI             1.1.3   2022-06-18 [1] CRAN (R 4.2.1)
#;-)  dbplyr          2.2.1   2022-06-27 [1] CRAN (R 4.2.1)
#;-)  digest          0.6.30  2022-10-18 [1] CRAN (R 4.2.1)
#;-)  dplyr         * 1.0.10  2022-09-01 [1] CRAN (R 4.2.1)
#;-)  ellipsis        0.3.2   2021-04-29 [1] CRAN (R 4.2.1)
#;-)  evaluate        0.18    2022-11-07 [1] CRAN (R 4.2.2)
#;-)  fansi           1.0.3   2022-03-24 [1] CRAN (R 4.2.1)
#;-)  fastmap         1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
#;-)  forcats       * 0.5.2   2022-08-19 [1] CRAN (R 4.2.1)
#;-)  fs              1.5.2   2021-12-08 [1] CRAN (R 4.2.0)
#;-)  gargle          1.2.1   2022-09-08 [1] CRAN (R 4.2.1)
#;-)  generics        0.1.3   2022-07-05 [1] CRAN (R 4.2.1)
#;-)  ggplot2       * 3.4.0   2022-11-04 [1] CRAN (R 4.2.2)
#;-)  glue            1.6.2   2022-02-24 [1] CRAN (R 4.2.1)
#;-)  googledrive     2.0.0   2021-07-08 [1] CRAN (R 4.2.0)
#;-)  googlesheets4   1.0.1   2022-08-13 [1] CRAN (R 4.2.1)
#;-)  gtable          0.3.1   2022-09-01 [1] CRAN (R 4.2.1)
#;-)  haven           2.5.1   2022-08-22 [1] CRAN (R 4.2.1)
#;-)  hms             1.1.2   2022-08-19 [1] CRAN (R 4.2.1)
#;-)  htmltools       0.5.4   2022-12-07 [1] CRAN (R 4.2.2)
#;-)  httr            1.4.4   2022-08-17 [1] CRAN (R 4.2.1)
#;-)  jsonlite        1.8.4   2022-12-06 [1] CRAN (R 4.2.2)
#;-)  knitr           1.41    2022-11-18 [1] CRAN (R 4.2.2)
#;-)  lifecycle       1.0.3   2022-10-07 [1] CRAN (R 4.2.1)
#;-)  lubridate       1.9.0   2022-11-06 [1] CRAN (R 4.2.2)
#;-)  magrittr        2.0.3   2022-03-30 [1] CRAN (R 4.2.1)
#;-)  modelr          0.1.10  2022-11-11 [1] CRAN (R 4.2.2)
#;-)  munsell         0.5.0   2018-06-12 [1] CRAN (R 4.2.1)
#;-)  pillar          1.8.1   2022-08-19 [1] CRAN (R 4.2.1)
#;-)  pkgconfig       2.0.3   2019-09-22 [1] CRAN (R 4.2.1)
#;-)  purrr         * 0.3.5   2022-10-06 [1] CRAN (R 4.2.1)
#;-)  R6              2.5.1   2021-08-19 [1] CRAN (R 4.2.1)
#;-)  Rcpp            1.0.9   2022-07-08 [1] CRAN (R 4.2.1)
#;-)  readr         * 2.1.3   2022-10-01 [1] CRAN (R 4.2.1)
#;-)  readxl          1.4.1   2022-08-17 [1] CRAN (R 4.2.1)
#;-)  reprex          2.0.2   2022-08-17 [1] CRAN (R 4.2.1)
#;-)  rlang           1.0.6   2022-09-24 [1] CRAN (R 4.2.1)
#;-)  rmarkdown       2.18    2022-11-09 [1] CRAN (R 4.2.2)
#;-)  rstudioapi      0.14    2022-08-22 [1] CRAN (R 4.2.1)
#;-)  rvest           1.0.3   2022-08-19 [1] CRAN (R 4.2.1)
#;-)  scales          1.2.1   2022-08-20 [1] CRAN (R 4.2.1)
#;-)  sessioninfo     1.2.2   2021-12-06 [1] CRAN (R 4.2.1)
#;-)  stringi         1.7.8   2022-07-11 [1] CRAN (R 4.2.2)
#;-)  stringr       * 1.5.0   2022-12-02 [1] CRAN (R 4.2.2)
#;-)  tibble        * 3.1.8   2022-07-22 [1] CRAN (R 4.2.1)
#;-)  tidyr         * 1.2.1   2022-09-08 [1] CRAN (R 4.2.1)
#;-)  tidyselect      1.2.0   2022-10-10 [1] CRAN (R 4.2.1)
#;-)  tidyverse     * 1.3.2   2022-07-18 [1] CRAN (R 4.2.2)
#;-)  timechange      0.1.1   2022-11-04 [1] CRAN (R 4.2.2)
#;-)  tzdb            0.3.0   2022-03-28 [1] CRAN (R 4.2.0)
#;-)  utf8            1.2.2   2021-07-24 [1] CRAN (R 4.2.1)
#;-)  vctrs           0.5.1   2022-11-16 [1] CRAN (R 4.2.2)
#;-)  withr           2.5.0   2022-03-03 [1] CRAN (R 4.2.1)
#;-)  xfun            0.34    2022-10-18 [1] CRAN (R 4.2.1)
#;-)  xml2            1.3.3   2021-11-30 [1] CRAN (R 4.2.0)
#;-)  yaml            2.3.6   2022-10-18 [1] CRAN (R 4.2.1)
#;-) 
#;-)  [1] /home/georges/R/x86_64-pc-linux-gnu-library/4.2
#;-)  [2] /usr/local/lib/R/site-library
#;-)  [3] /usr/lib/R/site-library
#;-)  [4] /usr/lib/R/library
#;-) 
#;-) ──────────────────────────────────────────────────────────────────────────────
```

</details>
