
c->e[0] = a->e[0]*b->e[0] - 
            s*t*(a->e[1]*b->e[1] + a->e[2]*b->e[2] + a->e[3]*b->e[3]);
c->e[1] = t*a->e[0]*b->e[1] + s*a->e[1]*b->e[0] + 
            s*t*(-a->e[2]*b->e[3] + a->e[3]*b->e[2]);
c->e[2] = t*a->e[0]*b->e[2] + s*t*a->e[1]*b->e[3] + 
            s*a->e[2]*b->e[0] - s*t*a->e[3]*b->e[1];
c->e[3] = t*a->e[0]*b->e[3] - s*t*a->e[1]*b->e[2] + 
            s*t*a->e[2]*b->e[1] + s*a->e[3]*b->e[0];

