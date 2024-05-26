typedef unsigned char uint8_t;
typedef float float32_t;
typedef double float64_t;

struct data{
 int32_t tid;
 float32_t *x, *profile;
};

typedef struct data data_t;

//
#define ARG_BYTE(av, i) (*((uint8_t *)((av)[i])))
#define ARG_SHORT(av, i) (*((short *)((av)[i])))
#define ARG_INT(av, i) (*((int32_t *)((av)[i])))
#define ARG_FLOAT(av, i) (*((float32_t *)((av)[i])))
#define ARG_FLOAT64(av, i) (*((float64_t *)((av)[i])))
//
#define ARG_BYTE_ARRAY(av, i) ((uint8_t *)((av)[i]))
#define ARG_SHORT_ARRAY(av, i) ((short *)((av)[i]))
#define ARG_INT_ARRAY(av, i) ((int32_t *)((av)[i]))
#define ARG_FLOAT_ARRAY(av, i) ((float32_t *)((av)[i]))
#define ARG_FLOAT64_ARRAY(av, i) ((float64_t *)((av)[i]))
