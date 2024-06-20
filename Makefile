md:
		Rscript -e "devtools::build_readme()"

site: md
		Rscript -e "pkgdown::build_site()"

check:
		Rscript -e "devtools::check()"

clean:
		Rscript -e "devtools::clean_dll()"

checkfast:
		Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

test:
		Rscript -e "devtools::test()"

doc:
		Rscript -e "devtools::document()"

build:
		Rscript -e "devtools::build()"

install:
		Rscript -e "devtools::install()"

builddata: install
		Rscript "data-raw/mcmc_results.R"

buildfast:
		Rscript -e "devtools::build(vignettes = FALSE)"

style:
		Rscript -e "styler::style_pkg()"

pr: clean style check site
