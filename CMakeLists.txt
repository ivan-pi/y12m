cmake_minimum_required(VERSION 3.12)

project(y12m
        VERSION 1.0.0
        LANGUAGES Fortran C)

add_library(y12m src/y12m.f90
    src/legacy/y12cck.f
    src/legacy/y12mae.f
    src/legacy/y12maf.f
    src/legacy/y12mbe.f
    src/legacy/y12mbf.f
    src/legacy/y12mce.f
    src/legacy/y12mcf.f
    src/legacy/y12mde.f
    src/legacy/y12mdf.f
    src/legacy/y12mfe.f
    src/legacy/y12mge.f
    src/legacy/y12mhe.f
    src/legacy/y12mhf.f)
