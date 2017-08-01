/* Some functions that are used within raw_processing.c
 * This file is directly #included in raw_processing.c */


/* Measurements taken from 5D Mark II RAW photos using EXIFtool, surely Canon can't be wrong about WB mutipliers? */
static const int wb_kelvin[]   = {  2500,  3000,  3506,  4000,  4503,  5011,  5517,  6018,  6509,  7040,  7528,  8056,  8534,  9032,  9531, 10000 };
static const double wb_red[]   = { 1.349, 1.596, 1.731, 1.806, 1.954, 2.081, 2.197, 2.291, 2.365, 2.444, 2.485, 2.528, 2.566, 2.612, 2.660, 2.702 };
static const double wb_green[] = { 1.137, 1.112, 1.056, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 };
static const double wb_blue[]  = { 3.985, 3.184, 2.524, 2.103, 1.903, 1.760, 1.641, 1.542, 1.476, 1.414, 1.390, 1.363, 1.333, 1.296, 1.263, 1.229 };

static const double id_matrix[] = {1,0,0,0,1,0,0,0,1};

static const double xyz_to_rgb[] = {
    3.240710, -0.9692580,  0.0556352,
   -1.537260,  1.8759900, -0.2039960,
   -0.498571,  0.0415557,  1.0570700
};

/* Cone space! */
static const double ciecam02[] = {
    0.7328,  0.4296, -0.1624,
   -0.7036,  1.6975,  0.0061,
    0.0030,  0.0136,  0.9834
};

/* Returns multipliers for white balance by (linearly) interpolating measured 
 * Canon values... stupidly simple, also range limited to 2500-10000 (pls obey) */
void get_kelvin_multipliers_rgb(double kelvin, double * multiplier_output)
{
    int k = 0;

    while (1 < 2)
    {
        if ( kelvin > wb_kelvin[k] && kelvin < wb_kelvin[k+1] ) 
        {
            break;
        }

        /* Doesn't need interpolation */
        else if ( kelvin == wb_kelvin[k] )
        {
            multiplier_output[0] = wb_red[k];
            multiplier_output[1] = wb_green[k];
            multiplier_output[2] = wb_blue[k];
            return;
        }

        k++;
    }

    double diff1 = (double)wb_kelvin[k+1] - (double)wb_kelvin[k];
    double diff2 = (double)wb_kelvin[k+1] - kelvin;

    /* Weight between the two */
    double weight1 = diff2 / diff1;
    double weight2 = 1.0 - weight1;

    multiplier_output[0] = (   wb_red[k] * weight1) + (   wb_red[k+1] * weight2);
    multiplier_output[1] = ( wb_green[k] * weight1) + ( wb_green[k+1] * weight2);
    multiplier_output[2] = (  wb_blue[k] * weight1) + (  wb_blue[k+1] * weight2);
}

void get_kelvin_multipliers_xyz(double kelvin, double * multiplier_output)
{
    double xyz_wb[3];
    double target[3];
    get_kelvin_multipliers_rgb(8000.0, xyz_wb);
    get_kelvin_multipliers_rgb(kelvin, target);
    for (int i = 0; i < 3; ++i) multiplier_output[i] = target[i] / xyz_wb[i];
    multiplier_output[1] *= 0.75;
}

/* now fix for XYZ :[ */
    // for (int i = 0; i < 3; ++i) multiplier_output[i] -= 1.0;
    // multiplier_output[0]
    // for (int i = 0; i < 3; ++i) multiplier_output[i] += 1.0;



/* Adds contrast(S-curve) to a value in a unique(or not) way */
double add_contrast ( double pixel, /* Range 0.0 - 1.0 */
                      double dark_contrast_range, 
                      double dark_contrast_factor, 
                      double light_contrast_range, 
                      double light_contrast_factor )
{
    /* Adjust based on the range, to make it not seem like range also controls strength */
    double dc_factor = dark_contrast_factor / (dark_contrast_range * 2);
    double lc_factor = light_contrast_factor / (light_contrast_range * 2);

    if (pixel < dark_contrast_range) 
        pixel = pow(pixel, 1 + (dark_contrast_range - pixel) 
        * (dc_factor
        /* Effect is reduced closer to end of range for a smooth curve: */
        * (1.0 - pixel/dark_contrast_range)) );

    pixel = 1.0 - pixel; /* Invert for contrasting the highlights */

    if (pixel < light_contrast_range && pixel > 0) 
        pixel = pow(pixel, 1 + (light_contrast_range - pixel) 
        * (lc_factor 
        * (1.0 - pixel/light_contrast_range)) );

    pixel = 1.0 - pixel; /* Invert back */

    return pixel;
}


/* Calculates the final matrix (processing->main_matrix), and precalculates all the values;
 * Combines: exposure + white balance + camera specific adjustment all in one matrix! */
void processing_update_matrices(processingObject_t * processing)
{
    /* Temporary working matrix */
    double temp_matrix_a[9];
    double temp_matrix_b[9];
    double temp_matrix_c[9];

    /* (for shorter code) */
    int32_t ** pm = processing->pre_calc_matrix;

    /* Create a camera to XYZ matrix in temp_matrix_a */
    invertMatrix(processing->xyz_to_cam_matrix, temp_matrix_a);

    /* whitebalance */

    /* Convert XYZ to eye cone space -> to temp_matrix_b */
    // multiplyMatrices(temp_matrix_a, (double *)ciecam02, temp_matrix_b); /* No ciecam for now */
    multiplyMatrices(temp_matrix_a, (double *)id_matrix, temp_matrix_b);
    memcpy(temp_matrix_a, temp_matrix_b, 9 * sizeof(double));

    /* Multiply channels */
    for (int i = 0; i < 3; ++i) temp_matrix_b[i] *= processing->wb_multipliers[0];
    for (int i = 3; i < 6; ++i) temp_matrix_b[i] *= processing->wb_multipliers[1];
    for (int i = 6; i < 9; ++i) temp_matrix_b[i] *= processing->wb_multipliers[2];

    /* Convert back to XYZ space from cone space -> to temp_matrix_a */
    invertMatrix((double *)ciecam02, temp_matrix_c);
    // multiplyMatrices(temp_matrix_b, temp_matrix_c, temp_matrix_a); /* No ciecam for now */
    multiplyMatrices(temp_matrix_b, (double *)id_matrix, temp_matrix_a);

    /* Multiply the currently XYZ matrix back to RGB in to final_matrix */
    multiplyMatrices( temp_matrix_a,
                      processing->xyz_to_rgb_matrix,
                      processing->final_matrix );

    /* Exposure, done here if smaller than 0, or no tonemapping - else done at gamma function */
    if (processing->exposure_stops < 0.0 || !processing->tone_mapping)
    {
        double exposure_factor = pow(2.0, processing->exposure_stops);
        for (int i = 0; i < 9; ++i) processing->final_matrix[i] *= exposure_factor;
    }

    /* Matrix stuff done I guess */

    /* Precalculate 0-65535 */
    for (int i = 0; i < 9; ++i)
    {
        for (int j = 0; j < 65536; ++j)
        {
            pm[i][j] = (int32_t)((double)j * processing->final_matrix[i]);
        }
    }

    /* Highest green value - pixels at this value will need to be reconstructed */
    processing->highest_green = processing->pre_calc_gamma[ LIMIT16(MAX(pm[3][65535],pm[3][0]) + MAX(pm[4][65535],pm[4][0]) + MAX(pm[5][65535],pm[5][0])) ];

    /* This is nice */
    printMatrix(processing->final_matrix);

    /* done? */
}