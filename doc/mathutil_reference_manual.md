# mathutil reference manual
## Constants

The values of the following constants vary by configuration.
```
real_max
real_dig
real_mant_dig
real_max_exp
real_min_exp
real_max_10_exp
real_min_10_exp
```
The values of the following constants is fixed.
```
#define r_pi (real_t)3.141592653589793238462643383279502884
#define r_1 (real_t)1
#define r_epsilon (real_t)0.000001
```

## Type definitions
### Scalar
```
typedef real_t *real_p;
```
### Vector
```
typedef struct vec4_struct
{
	real_t x, y, z, w;
}vec4_t, *vec4_p;
```
### Quaternion
```
typedef vec4_t quat_t, *quat_p;
```
### Matrix
```
typedef struct mat4_struct
{
	vec4_t x, y, z, w;
}mat4_t, *mat4_p;
```

## Functions
### [Scalar functions](#scalar_functions)
```
real_t r_rnd (uint32_t *p_seed);
real_t r_sin (real_t x);
real_t r_cos (real_t x);
real_t r_tan (real_t x);
real_t r_abs (real_t x);
real_t r_sgn (real_t x);
real_t r_floor (real_t x);
real_t r_ceil (real_t x);
real_t r_atan (real_t x);
real_t r_exp (real_t x);
real_t r_log (real_t x);
real_t r_pow(real_t x, real_t y);
real_t r_max(real_t x, real_t y);
real_t r_min(real_t x, real_t y);
real_t r_mod(real_t x, real_t y);
real_t r_atan2 (real_t x, real_t y);
real_t r_clamp (real_t n, real_t min, real_t max);
real_t r_lerp (real_t a, real_t b, real_t s);
real_t r_hermite (real_t s);
real_t r_slerp (real_t a, real_t b, real_t s);
```
### [Vector functions](#vector_functions)
```
vec4_t vec4(real_t x, real_t y, real_t z, real_t w);
vec4_t vec4_abs(vec4_t v);
vec4_t vec4_sgn(vec4_t v);
vec4_t vec4_invert(vec4_t v);
real_t vec4_length(vec4_t v);
vec4_t vec4_normalize(vec4_t v);
vec4_t vec4_scale(vec4_t v, real_t s);
vec4_t vec4_pow(vec4_t v);
real_t vec4_dot(vec4_t v1, vec4_t v2);
vec4_t vec4_add(vec4_t v1, vec4_t v2);
vec4_t vec4_sub(vec4_t v1, vec4_t v2);
vec4_t vec4_mul(vec4_t v1, vec4_t v2);
vec4_t vec4_div(vec4_t v1, vec4_t v2);
vec4_t vec4_min(vec4_t v1, vec4_t v2);
vec4_t vec4_max(vec4_t v1, vec4_t v2);
vec4_t vec4_cross3(vec4_t v1, vec4_t v2);
vec4_t vec4_mul_mat4(vec4_t v, mat4_t m);
vec4_t vec4_mul_mat4_transpose(vec4_t v, mat4_t m);
vec4_t vec4_lerp(vec4_t v1, vec4_t v2, real_t s);
vec4_t vec4_slerp(vec4_t v1, vec4_t v2, real_t s);
```
### [Quaternion functions](#quaternion_functions)
```
quat_t quat(real_t x, real_t y, real_t z, real_t w);
quat_t quat_rot_axis(vec4_t axis, real_t angle);
quat_t quat_mul(quat_t q1, quat_t q2);
quat_t quat_add_vec(quat_t q, vec4_t v, real_t s);
```
### [Matrix functions](#matrix_functions)
```
mat4_t mat4(vec4_t x, vec4_t y, vec4_t z, vec4_t w);
mat4_t mat4_rot_x(real_t angle);
mat4_t mat4_rot_y(real_t angle);
mat4_t mat4_rot_z(real_t angle);
mat4_t mat4_rot_axis(vec4_t axis, real_t angle);
mat4_t mat4_rot_euler(real_t yaw, real_t pitch, real_t roll);
mat4_t mat4_from_quat(quat_t q);
mat4_t mat4_from_quat_transpose(quat_t q);
mat4_t mat4_transpose(mat4_t m);
mat4_t mat4_add(mat4_t l, mat4_t r);
mat4_t mat4_add_s(mat4_t m, real_t s);
mat4_t mat4_add_transpose(mat4_t m, real_t s);
mat4_t mat4_sub(mat4_t l, mat4_t r);
mat4_t mat4_sub_s(mat4_t m, real_t s);
mat4_t mat4_sub_transpose(mat4_t m, real_t s);
mat4_t mat4_mul(mat4_t l, mat4_t r);
mat4_t mat4_mul_s(mat4_t m, real_t s);
mat4_t mat4_mul_transpose(mat4_t l, mat4_t r);
```

## Scalar functions <a name="scalar_functions"></a>
#### r_rnd <a name="r_rnd"></a>
Pseudo random number generator (aka. PRNG).

    real_t r_rnd (uint32_t *p_seed);
##### Parameters:
p_seed: A pointer to update a seed for random number generation.
##### Remarks:
The function works the same as MSVC *rand()* except that the seed is not a global variable, that allows the instancing purpose of using the PRNG.

---
#### r_sin <a name="r_sin"></a>
Sine function in trigonometric functions.

    real_t r_sin (real_t x);
##### Remarks:
The function is somewhat faster than MSVC's *sin()* if the alternative SIMD function is available, but lower accuracy.

---
#### r_cos <a name="r_cos"></a>
Cosine function in trigonometric functions.

    real_t r_cos (real_t x);
##### Remarks:
The function is somewhat faster than MSVC's *cos()* if the alternative SIMD function is available, but lower accuracy.

---
#### r_tan <a name="r_tan"></a>
Tangent function in trigonometric functions.

    real_t r_tan (real_t x);
##### Remarks:
The function is somewhat faster than MSVC's *tan()* if the alternative SIMD function is available, but lower accuracy.

---
#### r_abs <a name="r_abs"></a>
Absolute value function.

    real_t r_abs (real_t x);
##### Remarks:
Generic absolute value function. Returns -x if x is less than 0, otherwise returns x.

---
#### r_sgn <a name="r_sgn"></a>
Get the sign of the input number.

    real_t r_sgn (real_t x);
##### Remarks:
Returns -1 if x is less than 0; returns 1 if x is greater than 0. If x equals to 0, it returns 0.

---
#### r_floor <a name="r_floor"></a>
Returns the largest integer value not greater than x.

    real_t r_floor (real_t x);
##### Remarks:
For example, r_floor(3.5) returns 3, and r_floor(-3.5) returns -4, **instead of -3**.

---
#### r_ceil <a name="r_ceil"></a>
Returns the smallest integer value not less than x.

    real_t r_ceil (real_t x);
##### Remarks:
For example, r_ceil(3.5) returns 4, and r_ceil(-3.5) returns -3, **instead of -4**.

---
#### r_atan <a name="r_atan"></a>
Arctangent function in trigonometric functions.

    real_t r_atan (real_t x);

---
#### r_exp <a name="r_exp"></a>
Computes the e (Euler's number, 2.7182818) raised to the given power x.

    real_t r_exp (real_t x);

---
#### r_log <a name="r_log"></a>
Computes the natural (base e) logarithm of x.

    real_t r_log (real_t x);

---
#### r_pow <a name="r_pow"></a>
Computes the value of x raised to the power y.

    real_t r_pow(real_t x, real_t y)

---
#### r_max <a name="r_max"></a>
Get the maximum number of two numbers.

    real_t r_max(real_t x, real_t y)
##### Remarks:
The type of the parameters and the return value is *real_t*, which is generally not an integer type. That means, don't use this function to get the maximum number of two integers if you want the code faster. Because the conversion of a floating number to an integer can be **very slow**.

---
#### r_min <a name="r_min"></a>
Get the minimum number of two numbers.

    real_t r_min(real_t x, real_t y)
##### Remarks:
The type of the parameters and the return value is *real_t*, which is generally not an integer type. That means, don't use this function to get the minimum number of two integers if you want the code faster. Because the conversion of a floating number to an integer can be **very slow**.

---
#### r_mod <a name="r_mod"></a>
Computes the floating-point remainder of the division operation x/y.

    real_t r_mod(real_t x, real_t y)
##### Remarks:
The floating-point remainder of the division operation _x/y_ calculated by this function is exactly the value _x - n*y_, where *n* is x/y with its fractional part truncated.

The returned value has the same sign as _x_ and is less or equal to _y_ in magnitude.

---
#### r_atan2 <a name="r_atan2"></a>
Computes the arc tangent of y/x using the signs of arguments to determine the correct quadrant.

    real_t r_atan2 (real_t x, real_t y);
##### Remarks:
Returns the arc tangent of y/x in the range [-π, +π] radians.

---
#### r_clamp <a name="r_clamp"></a>
Constrain a value to lie between two further values

    real_t r_clamp (real_t n, real_t min, real_t max);
##### Remarks:
Returns the value of *x* constrained to the range *min* to *max*. The returned value is computed as *r_min(r_max(x, min), max)*.

---
#### r_lerp <a name="r_lerp"></a>
Performs a scalar linear interpolation calculation.

    real_t r_lerp (real_t a, real_t b, real_t s);
##### Remarks:
Returns _a + (b - a) * s_

---
#### r_hermite <a name="r_hermite"></a>
Hermite interpolation, which is smoother than linear interpolation.

    real_t r_hermite (real_t s);
##### Remarks:
The function modifies the *s* param of *r_lerp* for Hermite interpolation, can be used as the following example:
```
real_t value;
value = r_lerp(a, b, r_hermite(s));
```
In the example, *value* is Hermite interpolated between *a* and *b* by *s*.
About Hermite interpolation, see [Hermite interpolation](https://en.wikipedia.org/wiki/Hermite_interpolation)

---
#### r_slerp <a name="r_slerp"></a>
Performs a scalar Hermite interpolation calculation.

    real_t r_slerp (real_t a, real_t b, real_t s);
##### Remarks:
Returns _a + (b - a) * r_hermite(r_clamp(s, 0, 1))_.

---
### Vector functions <a name="vector_functions"></a>

---
#### vec4 <a name="vec4"></a>
The constructor of vec4_t. Returns a 4-d vector with specific components.

    vec4_t vec4(real_t x, real_t y, real_t z, real_t w);

---
#### vec4_abs <a name="vec4_abs"></a>
Get the absolute value of each component of vector *v*, and returns the result vector.

    vec4_t vec4_abs(vec4_t v);
##### Remarks:
Returns _vec4(r_abs(v.x), r_abs(v.y), r_abs(v.z), r_abs(v.w))_;

---
#### vec4_sgn <a name="vec4_sgn"></a>
Get the sign of each component of vector *v*, and returns the result vector. The sign is represented by -1, 1, or 0.

    vec4_t vec4_sgn(vec4_t v);
##### Remarks:
Returns _vec4(r_sgn(v.x), r_sgn(v.y), r_sgn(v.z), r_sgn(v.w))_;

---
#### vec4_invert <a name="vec4_invert"></a>
Returns a vector whose each component is the opposite of each component of vector *v*.

    vec4_t vec4_invert(vec4_t v);
##### Remarks:
Returns _vec4(-v.x, -v.y, -v.z, -v.w)_;

---
#### vec4_length <a name="vec4_length"></a>
Get the length of the vector *v*.

    real_t vec4_length(vec4_t v);
##### Remarks:
Returns _r_sqr(vec4_dot(v, v))_ where *vec4_dot* calculates the dot product of the vector *v*.

---
#### vec4_normalize <a name="vec4_normalize"></a>
Returns the normalized vector *v*. Normalize a vector means let the length of the vector equals to 1.

    vec4_t vec4_normalize(vec4_t v);
##### Remarks:
If the length of the vector *v* is less than r_epsilon (which generally equal to 0.000001), a zero-length vector is returned.

---
#### vec4_scale <a name="vec4_scale"></a>
Scales a vector.

    vec4_t vec4_scale(vec4_t v, real_t s);
##### Remarks:
Returns _vec4(v.x * s, v.y * s, v.z * s, v.w * s)_;

---
#### vec4_pow <a name="vec4_pow"></a>
Computes the value of each component of vector *v* raised to the power y, and returns the result vector.

    vec4_t vec4_pow(vec4_t v);
##### Remarks:
Returns _vec4(r_pow(v.x), r_pow(v.y), r_pow(v.z), r_pow(v.w))_;

---
#### vec4_dot <a name="vec4_dot"></a>
Determines the dot product of two vectors.

    real_t vec4_dot(vec4_t v1, vec4_t v2);
##### Remarks:
Returns _v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w_.

---
#### vec4_add <a name="vec4_add"></a>
Add two vectors.

    vec4_t vec4_add(vec4_t v1, vec4_t v2);

---
#### vec4_sub <a name="vec4_sub"></a>
Subtracts two vectors.

    vec4_t vec4_sub(vec4_t v1, vec4_t v2);

---
#### vec4_mul <a name="vec4_mul"></a>
Multiply each component of two vectors, and returns the result vector.

    vec4_t vec4_mul(vec4_t v1, vec4_t v2);

---
#### vec4_div <a name="vec4_div"></a>
Divides each component of two vectors, and returns the result vector.

    vec4_t vec4_div(vec4_t v1, vec4_t v2);

---
#### vec4_min <a name="vec4_min"></a>
Returns a vector that is made up of the smallest components of two vectors.

    vec4_t vec4_min(vec4_t v1, vec4_t v2);

---
#### vec4_max <a name="vec4_max"></a>
Returns a vector that is made up of the largest components of two vectors.

    vec4_t vec4_max(vec4_t v1, vec4_t v2);

---
#### vec4_cross3 <a name="vec4_cross3"></a>
Determines the cross-product of two **3D vectors**. Note that the *w* component of both vectors was ignored.

    vec4_t vec4_cross3(vec4_t v1, vec4_t v2);
##### Remarks:
The *w* component of the returned vector is zero.

---
#### vec4_mul_mat4 <a name="vec4_mul_mat4"></a>
Transforms a vector by a given matrix.

    vec4_t vec4_mul_mat4(vec4_t v, mat4_t m);

---
#### vec4_mul_mat4_transpose <a name="vec4_mul_mat4"></a>
Transforms a vector by a transposed matrix.

    vec4_t vec4_mul_mat4_transpose(vec4_t v, mat4_t m);

---
#### vec4_lerp <a name="vec4_lerp"></a>
Performs linear interpolation calculations between two vectors.

    vec4_t vec4_lerp(vec4_t v1, vec4_t v2, real_t s);
##### Remarks:
This function can also be used on quaternion interpolation. After doing interpolation on quaternion, it's better to normalize the result quaternion.

---
#### vec4_slerp <a name="vec4_slerp"></a>
Performs Hermite interpolation calculations between two vectors.

    vec4_t vec4_slerp(vec4_t v1, vec4_t v2, real_t s);
##### Remarks:
This function can also be used on quaternion interpolation. After doing interpolation on quaternion, it's better to normalize the result quaternion.

---
### Quaternion functions <a name="quaternion_functions"></a>

---
#### quat <a name="quat"></a>
The constructor of quat_t. Returns a 4-d vector with specific components.

    quat_t quat(real_t x, real_t y, real_t z, real_t w);

---
#### quat_rot_axis <a name="quat_rot_axis"></a>
Create a quaternion which represents the rotation around an axis.

    quat_t quat_rot_axis(vec4_t axis, real_t angle);

---
#### quat_mul <a name="quat_mul"></a>
Multiply two quaternions. 

    quat_t quat_mul(quat_t q1, quat_t q2);
##### Remarks:
The result represents the rotation *q1* followed by the rotation *q2* (_Out = q2 * q1_). The input quaternions must be normalized before calling this function, otherwise an incorrect result is returned. The multiplication of quaternions is not commutative.

To normalize a quaternion, use [vec4_normalize()](#vec4_normalize).

---
#### quat_add_vec <a name="quat_add_vec"></a>
Superimposes a rotation of a specific angle of rotation with the vector direction as an axis to a quaternion.

    quat_t quat_add_vec(quat_t q, vec4_t v, real_t s);
##### Remarks:
It's better to normalize the result quaternion.

---
### Matrix functions <a name="matrix_functions"></a>

---
#### mat4 <a name="mat4"></a>
The constructor of mat4_t, which represents a 4x4 matrix. The matrix have four 4-d vectors as it's content. The matrix is often used on vertex transformations.

    mat4_t mat4(vec4_t x, vec4_t y, vec4_t z, vec4_t w);

---
#### mat4_rot_x <a name="mat4_rot_x"></a>
Builds a matrix that represents the rotation around the x-axis.

    mat4_t mat4_rot_x(real_t angle);

---
#### mat4_rot_y <a name="mat4_rot_y"></a>
Builds a matrix that represents the rotation around the y-axis.

    mat4_t mat4_rot_y(real_t angle);

---
#### mat4_rot_z <a name="mat4_rot_z"></a>
Builds a matrix that represents the rotation around the z-axis.

    mat4_t mat4_rot_z(real_t angle);

---
#### mat4_rot_axis <a name="mat4_rot_axis"></a>
Builds a matrix that represents the rotation around the specific axis.

    mat4_t mat4_rot_axis(vec4_t axis, real_t angle);

---
#### mat4_rot_euler <a name="mat4_rot_euler"></a>
Builds a matrix with a specified yaw, pitch, and roll.

    mat4_t mat4_rot_euler(real_t yaw, real_t pitch, real_t roll);

---
#### mat4_from_quat <a name="mat4_from_quat"></a>
Builds a rotation matrix from a quaternion.

    mat4_t mat4_from_quat(quat_t q);
##### Remarks:
In order to get the proper result matrix, The input quaternion **must be normalized**.

---
#### mat4_from_quat_transpose <a name="mat4_from_quat_transpose"></a>
Builds a transposed rotation matrix from a quaternion.

    mat4_t mat4_from_quat_transpose(quat_t q);
##### Remarks:
The transposed rotation matrix can perform a transformation opposite to the quaternion rotation direction.
In order to get the proper result matrix, The input quaternion **must be normalized**.

---
#### mat4_transpose <a name="mat4_transpose"></a>
Transpose a matrix.

    mat4_t mat4_transpose(mat4_t m);

---
#### mat4_add <a name="mat4_add"></a>
Calculates the sum of each component of two matrices.

    mat4_t mat4_add(mat4_t l, mat4_t r);

---
#### mat4_add_s <a name="mat4_add_s"></a>
Adds a scalar to each component of the matrix and returns the result matrix.

    mat4_t mat4_add_s(mat4_t m, real_t s);

---
#### mat4_add_transpose <a name="mat4_add_transpose"></a>
Transposes the second matrix, then adds each component of the first matrix to the transposed matrix, then return the result matrix.

    mat4_t mat4_add_transpose(mat4_t m, real_t s);

---
#### mat4_sub <a name="mat4_sub"></a>
Subtracts each component of the first matrix to the second matrix and returns the result matrix.

    mat4_t mat4_sub(mat4_t l, mat4_t r);

---
#### mat4_sub_s <a name="mat4_sub_s"></a>
Subtracts a scalar to each component of the matrix and returns the result matrix.

    mat4_t mat4_sub_s(mat4_t m, real_t s);

---
#### mat4_sub_transpose <a name="mat4_sub_transpose"></a>
Transposes the second matrix, then subtracts each component of the first matrix to the transposed matrix, then return the result matrix.

    mat4_t mat4_sub_transpose(mat4_t m, real_t s);

---
#### mat4_mul <a name="mat4_mul"></a>
Matrix multiplication.

    mat4_t mat4_mul(mat4_t l, mat4_t r);

---
#### mat4_mul_s <a name="mat4_mul_s"></a>
Multiplies each component of a matrix to a scalar and returns the result matrix.

    mat4_t mat4_mul_s(mat4_t m, real_t s);

---
#### mat4_mul_transpose <a name="mat4_mul_transpose"></a>
Transposes the second matrix, then do matrix multiplication.

    mat4_t mat4_mul_transpose(mat4_t l, mat4_t r);

## Configurations

When compiling *mathutil*, the following prepocessor definitions changes the behavior of *mathutil* library.

    MATHUTIL_USE_DOUBLE
If it is defined to non-zero, real_t will be typedef to double. All related functions will treat real_t as double.
Otherwise, real_t will be typedef to float. By default, it is undefined and *mathutil* will typedef real_t to float.

    MATHUTIL_DETECT_SIMD
If defined to non-zero, call to any function of *mathutil* causes its initialization. *mathutil* will detect CPUID and use SIMD instructions which is compatible with the current CPU.
By default, *mathutil* detects the compiler. If the compiler is GCC or clang, *MATHUTIL_DETECT_SIMD* will be **defined to 0**. Otherwise *MATHUTIL_DETECT_SIMD* will be defined to 1.

    MATHUTIL_NOT_ALIGNED
    MATHUTIL_ASSUME_ALIGNED
The options above were ignored if *MATHUTIL_DETECT_SIMD* is 0. If any of them is defined to non-zero, any SIMD intrinsics that need memory address to be aligned (e.g. _mm_load_ps) will be used to replace the unaligned SIMD intrinsics.
The aligned intrinsics is somewhat faster than the unaligned one.