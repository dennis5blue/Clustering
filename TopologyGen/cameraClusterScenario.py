# To execute in the blender console
# filename = "/home/dennisyu/Documents/Clustering/TopologyGen/cameraClusterScenario.py"
# exec(compile(open(filename).read(), filename, 'exec'))

import bpy
import random
import math

###City settings
#city size 2000m^2, Seed = 4
#Create terrain = 1, Create meterials (sampling size) = 1, (texture size) = 2000, maximum height = 300, water level = 0, random options = 0
#Street texture size = 5000, slope limit = 5, Manhattan
#high poly models = 5, Total meterials = 15, texture size = 5000

# way of remove object, we need to select the objects first
# gather list of items of interest.
bpy.ops.object.select_all(action="DESELECT")
# remove all objects except the city.
removeList = [item.name for item in bpy.data.objects if item.type == "CAMERA" or item.type == "LAMP"]
for obj in removeList:
	bpy.data.objects[obj].select = True
bpy.ops.object.delete()

#move all the objects
bpy.ops.object.select_all(action="SELECT")
cityObjectList = [item.name for item in bpy.data.objects if item.select == True]
for obj in cityObjectList:
	bpy.data.objects[obj].location = [-2.5, -2.5, -3.0]

#bpy.ops.group.create(name="myCity")
#bpy.ops.object.group_link(group="myCity")
#myCity.location = [0.0, 0.0, 0.0]
#bpy.ops.object.location = [0.0, 0.0, 0.0]

# add light
scene = bpy.context.scene
lamp_data = bpy.data.lamps.new(name="sun1", type='SUN')
lamp_data.energy = 1.0
lamp_data.sky.use_sky = True
lamp_data.sky.use_atmosphere = True
lamp_object = bpy.data.objects.new(name="sun1", object_data=lamp_data)
scene.objects.link(lamp_object)
lamp_object.location = (0.0, 0.0, 10.0)
lamp_object.select = True

lamp_data = bpy.data.lamps.new(name="sun2", type='HEMI')
lamp_data.energy = 1.0
lamp_object = bpy.data.objects.new(name="sun2", object_data=lamp_data)
scene.objects.link(lamp_object)
lamp_object.location = (0.0, 0.0, 10.0)
lamp_object.select = True

scene.objects.active = lamp_object

#add camera (random position and sensing direction)
#city radius = 500m
for cam in range(20):
	cameraPos = (random.uniform(-2.5,2.5), random.uniform(-2.5,2.5), 0.1)
	cameraRotation = (math.pi*90.0/180, math.pi*0.0/180, -math.pi*random.randrange(0,360,10)/180)
	bpy.ops.object.camera_add(view_align=False, enter_editmode=False, location=cameraPos, rotation=cameraRotation);


bpy.data.scenes["Scene"].render.resolution_x = 1280*2;
bpy.data.scenes["Scene"].render.resolution_y = 720*2;
# object has a attribute named type which can query the type
# print(str(object.type))
for object in bpy.data.objects:
    print(object.name+" is at location " + str(object.location))

# to setup the camera field of view(FoV)
#bpy.data.scenes["Scene"].camera.data.angle_x
#bpy.data.scenes["Scene"].camera.data.angle_y
	
	
print('\nPrint Scenes...'); 
sceneKey = bpy.data.scenes.keys()[0]; 
print('Using Scene['  + sceneKey + ']');
c=0;
logFileName="/home/dennisyu/Documents/Clustering/sourceData/image_png/log.txt"
logFile=open(logFileName,"w+")
for obj in bpy.data.objects: 
# Find cameras that match cameraNames 
    if ( obj.type =='CAMERA'): 
        print("Rendering scene["+sceneKey+"] with Camera["+obj.name+"]")
        logFile.write("Rendering scene["+sceneKey+"] with Camera["+obj.name+"]\n") 
        location = obj.location
        rotation = obj.rotation_euler 
        logFile.write("Pos "+str(location.x)+" "+str(location.y)+" "+str(location.z)+" Dir "+str(rotation.x)+" "+str(rotation.y)+" "+str(rotation.z)+"\n")
        # Set Scenes camera and output filename 
        bpy.data.scenes[sceneKey].camera = obj 

        # Render Scene and store the scene 
        bpy.ops.render.render( write_still=True ); 
        RR = "Render Result";
        bpy.data.images[RR].save_render("/home/dennisyu/Documents/Clustering/sourceData/image_png/camera_"+str(c)+".png");
        c = c + 1; 

logFile.close()
