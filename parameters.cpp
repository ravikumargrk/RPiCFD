// Mesh dimensions
const double LenX = 0.090;
const double LenY = 0.060;

// Source dimensions
const double S1X_begin = 0.020;
const double S1X_end = 0.034;
const double S1Y_begin = 0.024;
const double S1Y_end = 0.038;

const double S2X_begin = 0.053;
const double S2X_end = 0.062;
const double S2Y_begin = 0.029;
const double S2Y_end = 0.038;

// Boundary values

const double he = 10;
const double Tinf = 300;

// Element dimensions (static mesh)
const size_t Nx = 18, Ny = 12;
const double DX = LenX/(double)Nx, DY = LenY/(double)Ny;

// Time frame parameters
const double max_time = 2;
const double DT = 1;
const double theta = 1; // weight of implicit method.