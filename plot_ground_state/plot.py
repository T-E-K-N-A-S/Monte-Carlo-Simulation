from PIL import Image
import os
import sys




csv_path = '/home/shikhar/Documents/sanket/Monte-Carlo-Simulation/csvs'
file_list =  os.listdir(csv_path)
file_list = [e for e in file_list if e.endswith('.csv') and e.startswith('Si')] 

for each_file in file_list:
    file_path = os.path.join(csv_path,each_file)
    f = open(file_path,'r')
    spin_mat =  f.readlines()

    im = Image.new("RGB",(128,128))
    pix = im.load()
    i = 0
    for row in spin_mat:
        j = 0
        for col in  row.strip().split(',')[:-1]:
            if col == '0' or col =='-1':
                pix[i,j] = (255,0,0)
                #data += chr(255) + chr(0) + chr(0) 
            else :
                
                pix[i,j] = (0,0,255)
                #data += chr(0) + chr(255)+ chr(0)
            #data += chr(255) + chr(255) + chr(255)
            j += 1
        i += 1
    #print i , j
    img_name = each_file  + ' plot.png'
    im.save(img_name,"PNG")
    print img_name
    f.close()
   
