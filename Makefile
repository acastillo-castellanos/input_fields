include $(BASILISK)/Makefile.defs

QCC = $(BASILISK)/qcc
CFLAGS += -O3 -DMTRACE=3 #-catch

all_tests: test_input_matrix.tst \
           test_init_fft.tst \
           test_init_fft2.tst \
           test_init_fft3.tst \
           test_init_fft4.tst \
           test_init_2Dto3D.tst \
					 test_init_2Dto3D.3D.tst

.PHONY: all_tests
