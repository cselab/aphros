
enum ISABELA_status
isabelaDeflateFloat(struct isabela_stream* ib_stream, enum ISABELA_flush ib_flush)
{
    uint32_t element_idx = 0;

    uint32_t element_byte_size = ib_stream->config->element_byte_size;
    uint32_t input_elements = ib_stream->avail_in / element_byte_size;

    uint32_t bits_per_index_element = 0;
    uint32_t ints_per_window = 0;
    uchar *in_buffer_uc = (uchar*) ib_stream->next_in;
    uchar *out_buffer_uc = (uchar*) ib_stream->next_out;

    uint32_t nwindows = 0;
    uint32_t window_idx = 0;

    uint32_t index_size = 0;
    uint32_t total_index_size = 0;
    size_t coefficients_size = 0;
    size_t error_size = 0;
    size_t comp_error_size = 0;
    size_t current_comp_error_size = 0;
    size_t current_decomp_error_size = 0;

    float *sorted_window_buffer = 0;
    float *approximate_window_buffer = 0;
    uint32_t *permuted_index = 0;
    void *error_tmp_buffer;

    if((ib_stream->avail_in % element_byte_size) != 0)
        return ISABELA_ERROR;

    print_isabela_config (ib_stream->config);
    // metadata is the following in binary format:
    //    4byte window_size, 4byte ncoefficients per window, 1byte data size, 4byte error rate
    // write metadata into the output stream first
    ib_stream->meta_ptr = out_buffer_uc;

    /*
     * Calculate offsets for index, coefficients and errors
     */
    nwindows                        = input_elements / ib_stream->config->window_size +
                                      (input_elements % ib_stream->config->window_size == 0? 0: 1);


    isabela_serialize_config (ib_stream->config, out_buffer_uc);
    out_buffer_uc += isabela_get_config_size (ib_stream->config);

    // npoints
    SET (uint32_t, out_buffer_uc, input_elements);
    out_buffer_uc += sizeof (uint32_t);

    // nwindows
    SET (uint32_t, out_buffer_uc, nwindows);
    out_buffer_uc                   += sizeof (uint32_t);

    bits_per_index_element          = calculate_bits_needed (ib_stream->config->window_size - 1);

    ints_per_window                 = (ib_stream->config->window_size * bits_per_index_element / (sizeof (uint32_t) * 8)) +
                                      ((ib_stream->config->window_size * bits_per_index_element % (sizeof (uint32_t) * 8)) == 0? 0: 1);

    ib_stream->error_size_ptr       = out_buffer_uc;

    // Advance the output buffer pointer
    // This is to store error and index offsets for each window
    out_buffer_uc += (2 * nwindows) * sizeof (uint32_t);

    // Perform ISABELA compression

    sorted_window_buffer        = (float *) malloc (sizeof (float) * ib_stream->config->window_size);
    approximate_window_buffer   = (float *) malloc (sizeof (float) * ib_stream->config->window_size);
    error_tmp_buffer            = (int32_t *) malloc (MAX_LEVEL * sizeof (int32_t) * ib_stream->config->window_size);

    for (window_idx = 0; window_idx < nwindows; window_idx ++) {

        uint32_t start               = (window_idx) * ib_stream->config->window_size;
        uint32_t end                 = MIN ((window_idx + 1) * ib_stream->config->window_size - 1, (input_elements - 1));
        uint32_t current_window_size = end - start + 1;

        error_size                = 0;
        index_size                = 0;

        memcpy (sorted_window_buffer,
                (float *) ib_stream->next_in + start,
                current_window_size * element_byte_size
                );

        // Sort based on window size, and build index
        sort_and_bit_index_float (sorted_window_buffer,
                                    element_byte_size,
                                    current_window_size,
                                    bits_per_index_element, 
                                    (uint32_t *) out_buffer_uc,
                                    &index_size
                                    );
        out_buffer_uc += index_size * sizeof (uint32_t);

        uint32_t compressed_buffer_size = 0;

        // B-Spline compression
        if (ib_stream->config->transform == ISABELA_BSPLINES) {
            compress_bspline_float   (sorted_window_buffer, current_window_size, (float *) out_buffer_uc, ib_stream->config->ncoefficients);
            compressed_buffer_size = ib_stream->config->ncoefficients * sizeof (float);
        } else if (ib_stream->config->transform == ISABELA_WAVELETS) {
            compress_wavelets_float (sorted_window_buffer, current_window_size, ib_stream->config->ncoefficients, out_buffer_uc, &compressed_buffer_size);
        }
        
        // Build the error coefficients
        if (ib_stream->config->error_rate != 0) {
            if (ib_stream->config->transform == ISABELA_BSPLINES) {
                decompress_bspline_float (approximate_window_buffer, current_window_size, (float *) out_buffer_uc, ib_stream->config->ncoefficients);
            } else {
                uint32_t wavelet_compressed_buffer_size = 0;
                decompress_wavelets_float (approximate_window_buffer, current_window_size, out_buffer_uc, &wavelet_compressed_buffer_size);
                assert (compressed_buffer_size == wavelet_compressed_buffer_size);
            }

            out_buffer_uc += compressed_buffer_size;

            for (element_idx = 0; element_idx < current_window_size; element_idx ++) {
                // Perform multi-level error encoding
                error_size        += error_encode_multilevel (sorted_window_buffer [element_idx],
                                                            approximate_window_buffer [element_idx],
                                                            (int32_t *) error_tmp_buffer + error_size,
                                                            ib_stream->config->error_rate
                                                            );
            }

            // Compress the error for each window,
            // and insert into the output buffer
            current_comp_error_size       = zlib2_compress_entire (error_size * sizeof (int32_t),
                                                        (char *) error_tmp_buffer,
                                                        error_size * sizeof (int32_t),
                                                        (char *) out_buffer_uc,
                                                        9
                                                        );
		    current_decomp_error_size = error_size * sizeof (int32_t);

            out_buffer_uc               += current_comp_error_size * sizeof (uchar);

        } else {
            out_buffer_uc += compressed_buffer_size;
            current_comp_error_size = 0;
            error_size = 0;
        }

        SET (uint32_t, (uint32_t *) ib_stream->error_size_ptr + 2 * window_idx, current_comp_error_size)
        SET (uint32_t, (uint32_t *) ib_stream->error_size_ptr + 2 * window_idx + 1, error_size * sizeof (int32_t))
    }

    // <window_size> <ncoefficient> <element size> <coefficient size> <error rate>
    // <nwindows>
    // <compressed error size 1> <compressed error size 2> ... <compressed error size NW> <-- Error size pointer
    // <index_1> <coefficients_1> <compressed_error_1>
    // <index_2> <coefficients_2> <compressed_error_2>
    // .
    // .
    // .
    // <index_NW> <coefficients_NW> <compressed_error_NW>

    free (sorted_window_buffer);
    free (approximate_window_buffer);
    free (error_tmp_buffer);

    ib_stream->end_ptr = (char *) out_buffer_uc;
    ib_stream->avail_out = (char *) ib_stream->end_ptr - (char *) ib_stream->meta_ptr;

    return ISABELA_SUCCESS;
}

enum ISABELA_status
isabelaInflateFloat(struct isabela_stream* ib_stream, enum ISABELA_flush ib_flush)
{
    uint32_t element_count;
    uint32_t element_idx;
    uint32_t i;

    uint32_t input_elements;
    uint32_t tmp_buffer_offset = 0;

    uchar *tmp_buf;
    uchar *tmp_ptr;

    uchar *in_buffer_uc = (uchar*) ib_stream->next_out;
    uchar *out_buffer_uc = (uchar*) ib_stream->next_in;

    uint32_t nwindows;
    uint32_t bits_per_index_element = 0;
    uint32_t ints_per_window = 0;

    uint32_t window_idx = 0;

    uint32_t index_size = 0;
    size_t coefficients_size = 0;
    size_t error_size = 0;
    size_t comp_error_size = 0;
    size_t current_comp_error_size = 0;
    size_t current_decomp_error_size = 0;

    float *original_window_buffer = 0;
    float *approximate_window_buffer = 0;
    uint32_t *permuted_index = 0;
    void *error_tmp_buffer;

    isabela_deserialize_config (ib_stream->config, out_buffer_uc);
    out_buffer_uc += isabela_get_config_size (ib_stream->config);

    // npoints
    GET (uint32_t, out_buffer_uc, input_elements);
    out_buffer_uc += sizeof (uint32_t);

    // nwindows
    GET (uint32_t, out_buffer_uc, nwindows);
    out_buffer_uc += sizeof (uint32_t);

    bits_per_index_element          = calculate_bits_needed (ib_stream->config->window_size - 1);

    ints_per_window                 = (ib_stream->config->window_size * bits_per_index_element / (sizeof (uint32_t) * 8)) +
                                      ((ib_stream->config->window_size * bits_per_index_element % (sizeof (uint32_t) * 8)) == 0? 0: 1);

    ib_stream->error_size_ptr       = (char *) out_buffer_uc;

    out_buffer_uc                   += (2 * nwindows) * sizeof (uint32_t);

    print_isabela_config (ib_stream->config);

    approximate_window_buffer   = (float *) malloc (sizeof (float) * ib_stream->config->window_size);
    permuted_index              = (uint32_t *) malloc (sizeof (uint32_t) * ib_stream->config->window_size);
    error_tmp_buffer            = (int32_t *) malloc (MAX_LEVEL * sizeof (int32_t) * ib_stream->config->window_size);

    // For each window
    for (window_idx = 0; window_idx < nwindows; window_idx ++) {

        uint32_t start               = (window_idx) * ib_stream->config->window_size;
        uint32_t end                 = MIN ((window_idx + 1) * ib_stream->config->window_size - 1, (input_elements - 1));
        uint32_t current_window_size = end - start + 1;

        tmp_buffer_offset       = 0;
        index_size                = 0;

        // uncompress index_w into permuted_index
        read_from_bitstream (current_window_size, bits_per_index_element, permuted_index, (uint32_t *) out_buffer_uc, &index_size);
        out_buffer_uc += index_size * sizeof (uint32_t);

        // <coefficients_1> <coefficients_2> ... <coefficients_NW> <-- Coefficients Pointer
        // B-Spline compression
        if (ib_stream->config->transform == ISABELA_BSPLINES) {
            decompress_bspline_float (approximate_window_buffer, current_window_size, (float *) out_buffer_uc, ib_stream->config->ncoefficients);
            out_buffer_uc           += ib_stream->config->ncoefficients * sizeof (float);
        } else if (ib_stream->config->transform == ISABELA_WAVELETS) {
            uint32_t wavelet_compressed_buffer_size = 0;
            decompress_wavelets_float (approximate_window_buffer, current_window_size, out_buffer_uc, &wavelet_compressed_buffer_size);
            out_buffer_uc += wavelet_compressed_buffer_size;
        }

        // <compressed error size 1> <compressed error size 2> ... <compressed error size NW> <-- Error size pointer
        GET (uint32_t, (uint32_t *) ib_stream->error_size_ptr + 2 * window_idx, current_comp_error_size);
        GET (uint32_t, (uint32_t *) ib_stream->error_size_ptr + 2 * window_idx + 1, current_decomp_error_size);

        if (current_comp_error_size > 0 && ib_stream->config->error_rate > 0) {
            // Compress the error for each window,
            // and insert into the output buffer
            zlib2_decompress_entire (current_comp_error_size,
                                    (char *) out_buffer_uc,
                                    current_decomp_error_size,
                                    (char *) error_tmp_buffer
                                    );
            out_buffer_uc += current_comp_error_size * sizeof (char);

            // Apply error into approximated
            for (i = 0; i < current_window_size; i ++) {
                approximate_window_buffer [i] = error_decode_multilevel (approximate_window_buffer[i], (int32_t *) error_tmp_buffer + tmp_buffer_offset, ib_stream->config->error_rate);
                tmp_buffer_offset += ((int32_t *) error_tmp_buffer) [tmp_buffer_offset] + 1;
            }
        }

        for (i = 0; i < current_window_size; i ++) {
            // Put data back into the original array
            ((float *) in_buffer_uc) [permuted_index [i]] = approximate_window_buffer [i];
        }

        in_buffer_uc = in_buffer_uc + sizeof (float) * current_window_size;
        ib_stream->avail_out += sizeof (float) * current_window_size;
    }

    ib_stream->end_ptr = out_buffer_uc;

    free (approximate_window_buffer);
    free (permuted_index);
    free (error_tmp_buffer);

    return ISABELA_SUCCESS;
}
