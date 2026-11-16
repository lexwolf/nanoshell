# Descriptive fits generated from corrected native-step pulse peak summaries.
# m is d(log10(Ipeak/Isat))/d(Gmax/Gth), fitted over 1.05--1.09.
omex_low_m = 99.533160676436879
omex_low_c = -7.6993154154085328
omex_low_alpha = 229.18357203214472
omex_high_a = 0.42378615442624457
omex_high_b = -0.42382883689529943
omex_low_fit(x) = (x>=1.05 ? 10**(omex_low_m*(x-1.05)+omex_low_c) : 1/0)
