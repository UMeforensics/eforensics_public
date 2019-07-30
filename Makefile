stepup: 
	Rscript test/set_up_test.r
helpfile1:
	Rscript test/helpfile_test.r
helpfile2:
	Rscript test/helpfile_correct_test.r
helpfile3:
	Rscript test/help_test.r
simulate:
	Rscript test/simulate_my_own_data_test.r
patrick:
	Rscript test/eforensics_test_20190126.R
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
dummy: test/efor_test1.Rmd
	Rscript -e 'rmarkdown::render("$<")'
clean:
	rm -rf *.html *.docx *.tar.gz figure/ cache/