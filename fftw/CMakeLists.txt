include_directories(${CMAKE_SOURCE_DIR}/fftw)

set(FFTW_C config.c executor.c fcr_1.c fcr_10.c fcr_11.c fcr_12.c fcr_128.c fcr_13.c fcr_14.c fcr_15.c fcr_16.c fcr_2.c fcr_3.c fcr_32.c fcr_4.c fcr_5.c fcr_6.c fcr_64.c fcr_7.c fcr_8.c fcr_9.c fftwf77.c fftwnd.c fhb_10.c fhb_16.c fhb_2.c fhb_3.c fhb_32.c fhb_4.c fhb_5.c fhb_6.c fhb_7.c fhb_8.c fhb_9.c fhf_10.c fhf_16.c fhf_2.c fhf_3.c fhf_32.c fhf_4.c fhf_5.c fhf_6.c fhf_7.c fhf_8.c fhf_9.c fn_1.c fn_10.c fn_11.c fn_12.c fn_13.c fn_14.c fn_15.c fn_16.c fn_2.c fn_3.c fn_32.c fn_4.c fn_5.c fn_6.c fn_64.c fn_7.c fn_8.c fn_9.c fni_1.c fni_10.c fni_11.c fni_12.c fni_13.c fni_14.c fni_15.c fni_16.c fni_2.c fni_3.c fni_32.c fni_4.c fni_5.c fni_6.c fni_64.c fni_7.c fni_8.c fni_9.c frc_1.c frc_10.c frc_11.c frc_12.c frc_128.c frc_13.c frc_14.c frc_15.c frc_16.c frc_2.c frc_3.c frc_32.c frc_4.c frc_5.c frc_6.c frc_64.c frc_7.c frc_8.c frc_9.c ftw_10.c ftw_16.c ftw_2.c ftw_3.c ftw_32.c ftw_4.c ftw_5.c ftw_6.c ftw_64.c ftw_7.c ftw_8.c ftw_9.c ftwi_10.c ftwi_16.c ftwi_2.c ftwi_3.c ftwi_32.c ftwi_4.c ftwi_5.c ftwi_6.c ftwi_64.c ftwi_7.c ftwi_8.c ftwi_9.c generic.c malloc.c planner.c putils.c rader.c rconfig.c rexec.c rexec2.c rfftwf77.c rfftwnd.c rgeneric.c rplanner.c timer.c twiddle.c wisdom.c wisdomio.c)

add_library(libfftw STATIC ${FFTW_C})


