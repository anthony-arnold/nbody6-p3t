
void vsub(const double l[3], const double r[3], double out[3]) {
    for (int i = 0; i < 3; i++) {
        out[i] = l[i] - r[i];
    }
}

void vsub_2d(const double l[2], const double r[2], double out[2]) {
    for (int i = 0; i < 2; i++) {
        out[i] = l[i] - r[i];
    }
}
