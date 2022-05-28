
void vsub(const double l[3], const double r[3], double out[3]) {
    for (int i = 0; i < 3; i++) {
        out[i] = l[i] - r[i];
    }
}
