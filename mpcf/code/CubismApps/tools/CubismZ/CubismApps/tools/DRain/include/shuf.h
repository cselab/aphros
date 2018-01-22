#ifdef __cplusplus
extern "C"
{
#endif
    size_t shuf_fw(void * const datastream,
			      const int nentries,
			      const int size_of_entry,
			      const int startbit, const int endbit);

    void shuf_bw(void * const datastream,
			     const int nentries,
			     const int size_of_entry,
			     const int startbit, const int endbit);
#ifdef __cplusplus
}
#endif
