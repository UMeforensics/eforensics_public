# Travis CI Testing
Everything in this folder is used to testing on Travis-ci.

## Adding more test cases

1. Add the test case into the /test folder, aka this directory
2. Inside the Makefile in the main directory folder, add the command that executes your test case, i.e.:

```
error_test:
	Rscript test/error_testing.R
```

where the first line is the makefile command to execute the terminal command in the second line,
thus to evoke this particle test case we would run:

```
make error_test
```

3. Add the makefile command to be ran into the .travis.yml file, under 'script:' section, i.e.:

```
script:
    ...
    travis_wait 90 make error_test
    ...
```

travis_wait command will force travis-ci to allow 90 mins for the test to complete.
There is no need to use travis_wait for smaller test cases, but for larger ones, tagging it in front is a safe choice

## Additional Notes
Inside the travis.yml file, I have added some comments to what each line does. In the 'before_install' section, some of the commands are encapsulated in if statements to seperate them out for each respective OS.

From what I've read, 'sudo: true' should be set almost always.

Caching packages are enabled too, it should be turned off periodically when running on the merge to master to ensure no compatibility issues.



## Limitations on Travis Tests
```
Runtime is longer than 120 mins
Limited to 3GB of Memory
Log output exceeds 4MB
```


```
Sources
https://github.com/softwaresaved/build_and_test_examples/blob/master/travis/README.md
```

## Description of Test Files
`vignette.R`: runs the code provided in the current help file

`help_files.R`: runs the help files and make sure they load properly

`error_test.R`: runs a series of datasets with common errors (missing abstension votes, missing total number of voters, missing covariates, etc.) to see if the package throws the proper errors

`toy_models.R`: runs a series of simulated datasets with the eforensics command

`king_county_pre`: uses King County precinct level data and processes the dataset such that it works with the eforensics command; demonstrates what fails and what needs to be corrected in order to properly run the function

`michigan_county`: uses Michigan County level data and processes the dataset such that it works with the eforensics command  
