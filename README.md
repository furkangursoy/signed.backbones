# Signed Backbone Extraction

This R package provides tools to extract the signed backbones of intrinsically dense weighted networks. 

Its counterpart in Python can be found at https://github.com/furkangursoy/signed_backbones.


## Installation

Install directly from this GitHub repository. Installation from CRAN will be available later.

```r
devtools::install_github("furkangursoy/signed.backbones")
```

## Example Usage

```r

karate_url <- 'https://raw.githubusercontent.com/furkangursoy/signed_backbones/main/examples/karate.txt'
karate_net <- read.csv(url(karate_url), header=FALSE, sep="\t")

karate_sbb <- signed.backbones::extract(karate_net, directed = FALSE, significance_threshold = 2.576, vigor_threshold = c(-0.1, 0.1))

```

## Citation

If you find this software useful in your work, please cite:

Furkan Gursoy and Bertan Badur. ["Extracting the signed backbone of intrinsically dense weighted networks"](https://arxiv.org/abs/2012.05216).



## Contributing

Please feel free to open an issue for bug reports, change requests, or other contributions.


## License

[MIT](https://choosealicense.com/licenses/mit/)
