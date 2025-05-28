import math

x1_array = [42.406, 80.298, 118.873, 158.120, 198.028, 238.579, 279.751, 321.518, 363.847, 406.700, 450.035, 493.803, 537.950, 582.420, 627.149, 672.072, 717.120, 762.224, 807.313]
x2_array = []
y2_array = []
x3_array = []
y3_array = []

for x1_k in x1_array:
  x2_k = - (1/2) * x1_k
  x2_array.append(x2_k)
  y2_k = - math.sqrt(3) * x2_k
  y2_array.append(y2_k)
  x3_k = x2_k
  x3_array.append(x2_k)
  y3_k = - y2_k
  y3_array.append(y3_k)

print(x2_array)
print(y2_array)
print(x3_array)
print(y3_array) 
