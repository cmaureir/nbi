

    def confirm_orbits(nmp,mps,mu):
        eccs = np.zeros(nmp)
        semis = np.zeros(nmp)
        for i in range(nmp):
            mp = mps[i]
            x = mp[1:4]
            v = mp[4:7]

            # r magnitude
            x_mag = np.linalg.norm(x)

            # Angular momentum
            # j = r x v
            j = np.cross(x, v)

            # Runge-Lenz-vector
            # e = { (v x j) / (G * m) }  - { r / |r| }
            e = np.cross(v, j)/(G * mu) - x/x_mag

            # Eccentricity
            ecc = np.linalg.norm(e)

            # Semi-major axis
            # a = ( j * j ) / (G * m * | 1 - ecc^2 | )
            a = j.dot(j)/(G* mu * np.fabs(1-ecc**2))

            semis[i] = a
            eccs[i] = ecc
        return eccs, semis

    def print_profile(m,r,fname):
        """Print the density profile

        Arguments:
        - `m`:
        - `r`:
        - `fname`:
        """
        #find out what dr is
        ofile = open(fname,'w')
        dr = r[1] - r[0]
        volfac = dr*4.0*np.pi
        for i in range(1,len(r)):
            radius = r[i]
            mass = m[i]
            dens = mass /( radius*radius*volfac)
            outline = "{0:.6e} {1:.6e}\n".format(radius,dens)
            ofile.write(outline)
        ofile.close()
