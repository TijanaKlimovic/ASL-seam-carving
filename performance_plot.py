from matplotlib import pyplot as plt 
import numpy as np

plt.style.use('fivethirtyeight')
fig = plt.figure()
ax = plt.gca()
ax.set_facecolor('#D3D3D3')
plt.title ('Intel Core i7 8550U @1.8 GHz, Compiler: gcc, OS: Windows 10 \n [F/C]', loc = 'left', fontsize = 14)
'Intel Core i7 8550U @1.8 GHz, Compiler: gcc, OS: Windows 10'

#images used on x axis 
same_image_diff_height_same_rem = []
diff_image_diff_width_same_rem = []
diff_image_same_size_same_rem = []

#performance on y axis
performance1 = []
performance2 = []
performance3 = []


plt.plot(same_image_diff_height_same_rem, performance1, marker ="o", label = "1", linewidth=1.5, color='#000000')
plt.plot(diff_image_diff_width_same_rem, performance2, marker ="o", label = "1", linewidth=1.5, color='#000000')
plt.plot(diff_image_same_size_same_rem, performance3, marker ="o", label = "1", linewidth=1.5, color='#000000')



plt.xlabel('Width size')
plt.xlabel('Height size')
plt.xlabel('Images')
#plt.ylabel("Performance [F/C]")


#plt.legend()

#plot for the 3rd exercise 
t = np.arange(0.01, 20.0, 0.01)

#plt.semilogx(dev_x3, dev_y, basex = 2, marker ="o", label = "1", linewidth=1.5, color='#000000')

#plt.xticks ([2 ** i for i in range (4,23 + 1)], fontsize = 13)
plt.xticks (same_image_diff_size, fontsize = 13)

plt.xlabel ("n", fontsize = 14) 
plt.yticks (fontsize = 14) 


ax.yaxis.grid(True, color = "#FFFFFF")
ax.xaxis.grid(True, color = "#D3D3D3")

#plt.tight_layout()

plt.show()