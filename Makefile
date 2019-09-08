error_test:
	Rscript test/error_testing.R
extreme_simulation:
	Rscript test/extreme_simulation.R
help_files:
	Rscript test/help_files.R
king_county_pre:
	Rscript test/king_county_precinctlevel_example.R
michigan_county:
	Rscript test/michigan_county_example.R
vignette:
	Rscript test/vignette.R
toy_models:
	Rscript test/toy_models.R
clean:
	rm -rf *.html *.docx *.tar.gz figure/ cache/
