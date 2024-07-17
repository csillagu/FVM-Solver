import numpy as np

class Assembly():
    def __init__(self, base_polygon, polygon_list) -> None:
        self.polygon_list = polygon_list
        self.base_polygon = base_polygon

    def plot(self, ax, colorMap = None, labels =False):
        if colorMap is None:
            self.base_polygon.plot(ax, "white", labels)
            for p in self.polygon_list:
                    p.plot(ax, labels=labels)
        else:
            polygons = [self.base_polygon]+self.polygon_list
            for p in polygons:
                if p.name not in colorMap:
                    p.plot(ax, labels=labels)
                else:
                    p.plot(ax, color=colorMap[p.name], labels=labels)
    

    def assemble(self):
        #Assembles the Assembly, assigns names to lines, points and polygons
        self._namePolygons()
        self._nameLines()

    def getPolygon(self, point):
        for p in self.polygon_list+[self.base_polygon]:
            if p.isPointInside(point):
                return p.name
        return ""

    def getPolygonNames(self):
        names =[p.name for p in [self.base_polygon]+self.polygon_list]
        return names
    
    def getLineNames(self):
        names = [l.name for p in [self.base_polygon]+self.polygon_list for l in p.lines]
        return names
    
    def getBaseLineNames(self):
        return [l.name for l in self.base_polygon.lines]

    def _namePolygons(self):
        self.base_polygon.setName("Polygon_0")
        for p in range(len(self.polygon_list)):
            self.polygon_list[p].setName("Polygon_"+str(p+1))
    
    def _nameLines(self):
        previous_names, previous_line_num = self.base_polygon.setLineNames(0, [])
        #previous_line_num = len(previous_names)
        for p in range(len(self.polygon_list)):
            newlines,previous_line_num = self.polygon_list[p].setLineNames(previous_line_num, previous_names) 
            previous_names=previous_names+newlines
    
    def getCoincidentLineName(self, line):
        for p in [self.base_polygon]+self.polygon_list:
            linename = p.getCoincidentLineName(line=line)
            if linename is not None:
                return linename
        
        return None
    def getLineByName(self, name):
        line= [l for p in [self.base_polygon]+self.polygon_list for l in p.lines if l.name == name]
        return line
    
    def getInnerPolygonLines(self):
        lines = [l for p in self.polygon_list for l in p.lines]
        return lines

class Polygon():
    def __init__(self, point_list) -> None:
        self.points = point_list
        self.lines = [Line(self.points[i], self.points[(i+1)%len(self.points)]) for i, p in enumerate(self.points)]
        self.name = ""

    def getPointsXY(self):
        points_y = [point.y for point in self.points]
        points_x = [point.x for point in self.points]
        return points_x, points_y
        
    def plot(self, ax, color="#c8ffc8", labels = True):
        points_x, points_y = self.getPointsXY()

        ax.fill(points_x, points_y, color=color)
        if labels:
            ax.text(sum(points_x)/len(points_x), sum(points_y)/len(points_y), self.name, ha='center', va='center', color="#5b745b", 
                bbox={'facecolor':'white','alpha':0.6,'edgecolor':'none','pad':1})
        for l in self.lines:
            l.plot(ax,"k-", labels = labels)
    
    def setName(self, name):
        self.name = name
        
    def lineAlreadyExists(line, lines):
        if len(lines) == 0:
            return False, ""
        for lin in lines:
            if lin == line:
                return True, lin.name
        return False, ""
    
    def setLineNames(self,start_index, previousLines):
        next_name = start_index
        for p in range(len(self.lines)):
            ex, name = Polygon.lineAlreadyExists(self.lines[p], previousLines)
            if ex:
                self.lines[p].setName(name)
            else:
                self.lines[p].setName("Line_"+str(next_name))
                next_name+=1
        return self.lines, next_name
    
    def centerMass(self):
        points_x, points_y = self.getPointsXY()
        return Vector(sum(points_x)/len(points_x), sum(points_y)/len(points_y))

    def isPointInside(self, point):

        # source: https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/amp/
        polygon = self.points

        num_vertices = len(polygon)
        
        x, y = point.x, point.y

        inside = False
    

        # Store the first point in the polygon and initialize the second point

        p1 = polygon[0]
    

        # Loop through each edge in the polygon

        for i in range(1, num_vertices + 1):

            # Get the next point in the polygon

            p2 = polygon[i % num_vertices]
    

            # Check if the point is above the minimum y coordinate of the edge

            if y > min(p1.y, p2.y):

                # Check if the point is below the maximum y coordinate of the edge

                if y <= max(p1.y, p2.y):

                    # Check if the point is to the left of the maximum x coordinate of the edge

                    if x <= max(p1.x, p2.x):

                        # Calculate the x-intersection of the line connecting the point to the edge

                        x_intersection = (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x
    

                        # Check if the point is on the same line as the edge or to the left of the x-intersection

                        if p1.x == p2.x or x <= x_intersection:

                            # Flip the inside flag

                            inside = not inside
    

            # Store the current point as the first point for the next iteration

            p1 = p2
    

        # Return the value of the inside flag

        return inside
    
    def getCoincidentLineName(self, line):
        for l in self.lines:
            if l.isLineOnThisLine(line=line):
                return l.name
        
        return None



class Line():
    def __init__(self, p1, p2) -> None:
        assert(p1!=p2)
        self.p1 = p1
        self.p2 = p2
        self.name = ""
    
    def plot(self, ax, fmt="k-", labels = True):
        ax.plot([self.p1.x, self.p2.x], [self.p1.y, self.p2.y], fmt)
        if labels:
            ax.text((self.p1.x+self.p2.x)/2, (self.p1.y+self.p2.y)/2, self.name, ha='center', va='center', 
                bbox={'facecolor':'white','alpha':0.6,'edgecolor':'none','pad':1},color = "black")
    
    def setName(self, name):
        self.name=name
    
    def __eq__(self, other: object) -> bool:
        
        return (self.p1 == other.p1) and (self.p2 == other.p2)
    
    def isPointOnLine(self, p):
        a = self.p1
        b = self.p2
        crossproduct = (p.y - a.y) * (b.x - a.x) - (p.x - a.x) * (b.y - a.y)

        # compare versus epsilon for floating point values, or != 0 if using integers
        if abs(crossproduct) > 1e-3:
            return False

        dotproduct = (p.x - a.x) * (b.x - a.x) + (p.y - a.y)*(b.y - a.y)
        if dotproduct < -1e-3:
            return False

        squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y)
        if dotproduct > squaredlengthba:
            return False

        return True
    
    def isLineOnThisLine(self, line):
        #Checks if "line" is on this line
        return self.isPointOnLine(line.p1) and self.isPointOnLine(line.p2)
    
    def isPerpendicular(self, line2):
        vec1 = [self.p2.x-self.p1.x, self.p2.y-self.p1.y]
        vec2 = [line2.p2.x-line2.p1.x, line2.p2.y-line2.p1.y]

        return vec1[0]*vec2[0]+vec1[1]*vec2[1] == 0
    
    def getNormal(self):
        dx = self.p2.x-self.p1.x
        dy = self.p2.y-self.p1.y
        return Vector(dy, -dx)
    
    def getLength(self):
        dx = self.p2.x-self.p1.x
        dy = self.p2.y-self.p1.y
        return np.sqrt(dx**2+dy**2)
    def getDirectionVector(self):
        dx = self.p2.x-self.p1.x
        dy = self.p2.y-self.p1.y
        return Vector(dx, dy)
    

    
    def getCenter(self):
        return Vector(sum([self.p1.x, self.p2.x])/2, sum([self.p1.y, self.p2.y])/2)



class Vector():
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def getPosition(self):
        return (self.x, self.y)
    
    def plot(self, ax, fmt = "bx"):
        ax.plot(self.x, self.y,fmt)

    def plotAsVector(self, ax, vect_0 = (0,0), color = "b", scale = 1):
        ax.arrow(vect_0.x,vect_0.y,self.x*scale,self.y*scale, color=color)
    

    def __eq__(self, other: object) -> bool:
        return (self.x == other.x) and (self.y == other.y)
    def __sub__(self, other):
        #self-other
        return Vector(self.x-other.x, self.y-other.y)
    
    def __mul__(self, other):
        return self.x*other.x+self.y*other.y
    
    def toNpArray(self):
        return np.array([[self.x], [self.y]])

