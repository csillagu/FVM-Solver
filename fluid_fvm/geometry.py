class Assembly():
    def __init__(self, base_polygon, polygon_list) -> None:
        self.polygon_list = polygon_list
        self.base_polygon = base_polygon

    def plot(self, ax):
        self.base_polygon.plot(ax, "white")
        for p in self.polygon_list:
            p.plot(ax)
    

    def assemble(self):
    #Assembles the Assembly, assigns names to lines, points and polygons
        self._namePolygons()
        self._nameLines()


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

class Polygon():
    def __init__(self, point_list) -> None:
        self.points = point_list
        self.lines = [Line(self.points[i], self.points[(i+1)%len(self.points)]) for i, p in enumerate(self.points)]

    def getPointsXY(self):
        points_y = [point.y for point in self.points]
        points_x = [point.x for point in self.points]
        return points_x, points_y
        
    def plot(self, ax, color="#c8ffc8"):
        points_x, points_y = self.getPointsXY()

        ax.fill(points_x, points_y, color = color, edgecolor="black")
    
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



class Line():
    def __init__(self, p1, p2) -> None:
        assert(p1!=p2)
        self.p1 = p1
        self.p2 = p2
        self.name = ""
    
    def plot(self, ax, fmt="b-"):
        ax.plot([self.p1.x, self.p2.x], [self.p1.y, self.p2.y])
        ax.text((self.p1.x+self.p2.x)/2, (self.p1.y+self.p2.y)/2, self.name)
    
    def setName(self, name):
        self.name=name
    
    def __eq__(self, other: object) -> bool:
        
        return (self.p1 == other.p1) and (self.p2 == other.p2)


class Point():
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def getPosition(self):
        return (self.x, self.y)
    
    def plot(self, ax, fmt = "bx"):
        ax.plot(self.x, self.y,fmt)

    def __eq__(self, other: object) -> bool:
        return (self.x == other.x) and (self.y == other.y)

