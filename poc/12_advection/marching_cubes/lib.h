struct Vec {
    float x;
    float y;
    float z;
};
struct Vec;
float GetOffset(float, float, float);
void GetColor(struct Vec*, struct Vec*);
void PrintHelp(void);
void GetNormal(struct Vec*, float, float, float, double(*f)(double, double, double, void*), void*);
	  
