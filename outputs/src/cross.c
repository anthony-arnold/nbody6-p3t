
void cross(const double l[3], const double r[3], double out[3]) {
    out[0] = l[1]*r[2] - l[2]*r[1];
    out[1] = l[2]*r[0] - l[0]*r[2];
    out[2] = l[0]*r[1] - l[1]*r[0];
}
