from math import *
# do a fly(samth, -2,2,41, zpx, -4,4, 41,0.02) scan
# dzpx is totally zpx changes from negetive to positive
# dsamth is totally samth changes from negetive to positive

def center_of_rotation(dzpx, dsamth, samth):
    del_samx = dzpx / radians(dsamth) * sin(radians(samth))
    del_samz = - del_samx / tan(radians(samth))
    return del_samx, del_samz

# Example input values (adjust as needed)
dzpx = 2.612
dsamth = 3.28
samth = 20.4

del_samx, del_samz = center_of_rotation(dzpx, dsamth, samth)
print("samx should move:", del_samx)
print("samz should move:", del_samz)
print("footprint ratio is:", sin(radians(samth)))
