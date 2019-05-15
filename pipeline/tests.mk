
#####################################################################################################
# unit tests
#####################################################################################################

test_basic:
	@$(MAKE) pp_basic c=../config/template/basic.cfg BASE_OUTDIR=`pwd`/out/basic BASE_TMPDIR=`pwd`/out/basic_tmp HPIPE_DIR=`pwd`/..
test_skip_assembly:
	@$(MAKE) pp_basic c=../config/template/skip_assembly.cfg BASE_OUTDIR=`pwd`/out/skip_assembly BASE_TMPDIR=`pwd`/out/skip_assembly_tmp HPIPE_DIR=`pwd`/.
test_multi_sites:
	@$(MAKE) pp_basic c=../config/template/multi_sites.cfg BASE_OUTDIR=`pwd`/out/multi_sites BASE_TMPDIR=`pwd`/out/multi_sites_tmp HPIPE_DIR=`pwd`/..
all_tests: test_basic test_skip_assembly test_multi_sites
