
using System.Collections.Generic;
using UnityEngine;

public class PhysicsRect : RigidBody
{
    // Possible combinations of forces
    public static int NO_FORCE = 0;
    public static int TURNING_LEFT_FORCE = 1 << 0;
    public static int TURNING_RIGHT_FORCE = 1 << 1;
    public static int FORWARD_FORCE = 1 << 2;
    public static int BACKWARD_FORCE = 1 << 3;
    
    // Force magnitudes
    public static float LEFT_RIGHT_FORCE_MAG = 5000;
    public static float BACK_FORWARD_FORCE_MAG = 10000;
    
    Vector2 COM;
    Vector2 radial;
    Vector2 forcePosition;
    Vector2 velocityTotal;
    float momentOfInertia;
    float width;
    float height;
    Color colour;
    Vector2 force;
    float linearAccel;
    float angularAccel;
    public float dragCoefficient = 250;
    public List<PhysicsRect> childRects;
    private int currentForces;
    private double circleRadius;
    private Vector2[] initialVerts;
    private List<Pair<Vector2, Vector2>> edges;

    public PhysicsRect(float x, float y, float width, float height, Color colour, float mass, float rotation)
    {
        base(mass, new Vector2(x, y), rotation, new Vector2(0, 0), 0f);
        
        velocityTotal =  new Vector2(0,0);
        this.width = width;
        this.height = height;
        this.colour = colour;
        childRects = new List<PhysicsRect>();
        COM = new Vector2();
        force = new Vector2();
        radial = new Vector2();
        linearAccel = 0;
        momentOfInertia = 0;
        forcePosition = new Vector2();
                
        circleRadius = Mathf.Sqrt(Mathf.Pow(width / 2, 2) + Mathf.Pow(height / 2, 2));
        
        initialVerts = new Vector2[4];
        initialVerts[0] = new Vector2(-width / 2, height / 2);
        initialVerts[1] = new Vector2(-width / 2, -height / 2);
        initialVerts[2] = new Vector2(width / 2, -height / 2);
        initialVerts[3] = new Vector2(width / 2, height / 2);
        
        vertices = new Vector2[4];
        edges = new List<Pair<Vector2, Vector2>>(initialVerts.Length);
        for(int i = 0; i < vertices.length; i++)
        {
            vertices[i] = initialVerts[i];
        }
        
        for(int i = 0; i < vertices.length; i++)
        {
            Vector2 nextV = i == vertices.length - 1 ? vertices[0] : vertices[i + 1];
            Pair<Vector2, Vector2> edge = new Pair<Vector2, Vector2>(vertices[i], nextV);
            edges.Add(edge);
        }
        updateVertices();
    }

    public void reset(float x, float y, float width, float height, Color colour, float mass, float rotation)
    {
        position.x = x;
        position.y = y;
        velocity.x = 0.0f;
        velocity.y = 0.0f;
        angularVelocity = 0.0f;
        this.width = width;
        this.height = height;
        this.colour = colour;
        this.mass = mass;
        this.rotation = rotation;
        force = new Vector2();
        radial = new Vector2();
        linearAccel = 0;
        momentOfInertia = 0;
        forcePosition = Vector2.zero ;
        velocityTotal = new Vector2(0,0);
    }

    public Vector2 backForwardForceLocation()
    {
        /*Vector2 pos = new Vector2(position.x - width/2, position.y);
        pos -= position;
        pos.rotate(rotation);
        pos += position;*/

        //TODO: MainThruster
        return Vector2.zero;//pos;
    }

    public Vector2 turnLeftForceLocation()
    {
        /*Vector2 pos = new Vector2( position.x + width/2, position.y - height/2);
        pos.sub(position);
        pos.rotate(rotation);
        pos.add(position);*/

        //TODO: LeftTurnThruster:
        return Vector2.zero;//pos
    }

    public Vector2 turnRightForceLocation()
    {
        /*Vector2 pos = new Vector2( position.x + width/2, position.y + height/2);
        pos.sub(position);
        pos.rotate(rotation);
        pos.add(position);*/

        //TODO: RightTurnThruster
        return Vector2.zero;//pos;
    }

    public void addChild(PhysicsRect rect)
    {
        childRects.Add(rect);
    }

    public float getTotalMass()
    {
        float totalMass = this.mass;

        foreach (PhysicsRect child in childRects)
        {
            totalMass += child.mass;
        }

        return totalMass;
    }

    public void centralForceOn(bool forward)
    {
        currentForces |= forward ? FORWARD_FORCE: BACKWARD_FORCE;
    }

    public void centralForceOff(bool forward)
    {
        if (forward)
        {
            if ((currentForces & FORWARD_FORCE) != 0)
                currentForces ^= FORWARD_FORCE;
        }
        else
        {
            if ((currentForces & BACKWARD_FORCE) != 0)
                currentForces ^= BACKWARD_FORCE;
        }
    }

    public bool isForceOn(int force)
    {
        return ((currentForces & force) != 0);
    }

    public void rightForceOn()
    {
        currentForces |= TURNING_RIGHT_FORCE;
    }

    public void rightForceOff()
    {
        if ((currentForces & TURNING_RIGHT_FORCE) != 0)
            currentForces ^= TURNING_RIGHT_FORCE;
    }

    public bool isForceOn()
    {
        return currentForces != NO_FORCE;
    }

    public void leftForceOn()
    {
        currentForces |= TURNING_LEFT_FORCE;
    }

    public void leftForceOff()
    {
        if ((currentForces & TURNING_LEFT_FORCE) != 0)
            currentForces ^= TURNING_LEFT_FORCE;
    }

    public override float getMomentOfInertia()
    {
        float moi = 0;
        float Icm = (float)((mass * ((float)Mathf.Pow(width, 2) + (float)Mathf.Pow(height, 2))) / 12.0);
        float H = (float)(Mathf.Sqrt(Mathf.Pow(position.x - COM.x,2) + Mathf.Pow(position.y - COM.y,2)));
        float I = Icm + (float)(mass * Mathf.Pow(H, 2));

        moi += I;

        foreach (PhysicsRect rect in childRects)
        {
            Icm = (float)((rect.mass * ((float)Mathf.Pow(rect.width, 2) + (float)Mathf.Pow(rect.height, 2))) / 12.0);
            H = (float)(Mathf.Sqrt(Mathf.Pow(rect.position.x - COM.x,2) + Mathf.Pow(rect.position.y - COM.y,2)));
            I = Icm + (float)(rect.mass * Mathf.Pow(H, 2));

            moi += I;
        }
        return moi;
    }

    private float determineAngularAcceleration(float I)
    {
        Vector2 r = radial;
        Vector2 f = force;
        float t = Vector3.Cross(r, f).z;

        float angularAccel = 0.0f;
        angularAccel = t / I;

        return angularAccel;
    }

    public override double getBoundingCircleRadius()
    {
        return circleRadius;
    }
    
    public override Vector2[] getVertices()
    {
        return vertices;
    }
    
    public override List<Pair<Vector2, Vector2>> getEdges()
    {
        return edges;
    }

    /**
     * Updates the physics simulation.
     * 
     * NOTE: I'm still not entirely sure the order in which things should be updated.
     * @param time The amount of time (in seconds) that has elapsed since the last
     * update.
     */
    public override void update(float time)
    {
        forcePosition.x = 0;
        forcePosition.y = 0;
        int numForces = 0;
        force.y = 0;
        force.x = 0;

        // Determine which forces are active and average their positions and
        // magnitudes
        if ((currentForces & TURNING_RIGHT_FORCE) != 0)
        {
            numForces++;

            forcePosition.x += position.x + width/2;
            forcePosition.y += position.y + height/2;
            force.x += LEFT_RIGHT_FORCE_MAG;
        }

        if ((currentForces & TURNING_LEFT_FORCE) != 0)
        {
            numForces++;

            forcePosition.x += position.x + width/2;
            forcePosition.y += position.y - height/2;
            force.x += LEFT_RIGHT_FORCE_MAG;
        }

        if ((currentForces & FORWARD_FORCE) != 0)
        {
            numForces++;

            forcePosition.x += position.x - width/2;
            forcePosition.y += position.y;
            force.x += BACK_FORWARD_FORCE_MAG;
        }

        if((currentForces & BACKWARD_FORCE) != 0)
        {
            numForces++;

            forcePosition.x += position.x - width/2;
            forcePosition.y += position.y;
            force.x -= BACK_FORWARD_FORCE_MAG;
        }

        // Average the forces
        numForces = numForces == 0 ? 1 : numForces; // Make this 1 to avoid dividing by 0
        forcePosition.x /= numForces;
        forcePosition.y /= numForces;

        forcePosition.sub(position); //TODO: Here
        forcePosition.rotate(rotation);
        forcePosition.add(position);
        force.rotate(rotation);

        // Calculate COM
        float totalMass = mass;
        COM.x = position.x * mass;
        COM.y = position.y * mass;

        foreach (PhysicsRect rect in childRects)
        {
            COM.x += rect.position.x * rect.mass;
            COM.y += rect.position.y * rect.mass;
            totalMass += rect.mass;
        }
        if (totalMass > 0)
        {
            COM.x /= totalMass;
            COM.y /= totalMass;
        }

        // update vel, accel, and position
        radial.x = forcePosition.x - COM.x;
        radial.y = forcePosition.y - COM.y;

        Vector2 acceleration = new Vector2(force.x/totalMass, force.y/totalMass);
        this.momentOfInertia = getMomentOfInertia();

        float angularAcceleration = determineAngularAcceleration(momentOfInertia);

        // Determine Angular Velocity
        float angularForce = (angularAcceleration * totalMass);
        angularVelocity = 1/dragCoefficient * (angularForce - (float)Mathf.Pow(Mathf.E, -dragCoefficient * time/totalMass) * (angularForce - dragCoefficient * angularVelocity));

        float rotationThisFrame =  angularVelocity  * time * (float)(180/Mathf.PI);
        rotation += rotationThisFrame;

        // Determine Velocity
        linearAccel = acceleration.x;
        velocity.x += acceleration.x * time;
        velocity.y += acceleration.y * time;

        // update position
        COM.x += velocity.x * time;
        COM.y += velocity.y * time;

        position.x += velocity.x * time;
        position.y += velocity.y * time;
        position.sub(COM);
        position.rotate(rotationThisFrame);
        position.add(COM);

        foreach (PhysicsRect rect in childRects)
        {
            rect.position.x += velocity.x * time;
            rect.position.y += velocity.y * time;
            rect.rotation = rotation;
            rect.position.sub(COM);
            rect.position.rotate(rotationThisFrame);
            rect.position.add(COM);
        }

        forcePosition.x  += velocity.x * time;
        forcePosition.y  += velocity.x * time;
        forcePosition.sub(COM); //TODO: Here
        forcePosition.rotate(rotationThisFrame);
        forcePosition.add(COM);

        final float eToCoeff = (float)Mathf.Pow(Mathf.E, -dragCoefficient * time/totalMass);
        velocity.x = 1/dragCoefficient * (force.x - eToCoeff * (force.x - dragCoefficient * velocity.x));
        velocity.y = 1/dragCoefficient * (force.y - eToCoeff * (force.y - dragCoefficient * velocity.y));

        updateVertices();
    }

   public override void setAngularVelocity(float vel)
    {
//        float momentOfInertia = getMomentOfInertia();
//        this.momentOfInertia = momentOfInertia;
//
//        float angularAcceleration = determineAngularAcceleration(momentOfInertia);
//
//        // Determine Angular Velocity
//        float angularForce = (angularAcceleration * totalMass);
//
//        angularVelocity = 1/dragCoefficient * (angularForce - (float)Mathf.Pow(Mathf.E, -dragCoefficient * time/totalMass) * (angularForce - dragCoefficient * angularVelocity));
    }
    
    /**
     * Updates the box's vertex positions.
     */
    public override void updateVertices()    
    {
        for(int i = 0; i < vertices.Length; i++)
        {
            vertices[i].x = initialVerts[i].x;
            vertices[i].y = initialVerts[i].y;
            vertices[i].rotate(rotation);
            vertices[i].x += position.x;
            vertices[i].y += position.y;
        }
    }
    
    public override Vector2 getCentreOfMass()
    {
        return COM;
    }
}
