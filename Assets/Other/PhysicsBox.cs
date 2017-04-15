

using System.Collections.Generic;
using UnityEngine;
/**
*
* @author Shane
*/
public class PhysicsBox : RigidBody {
    
    // Collision information
    private double circleRadius;
    
    // Whether the box is elastic
    public bool isElastic;
    
    // The length of one side of the box
    public float sideLen;
    
    // The positions of each vertex relative the centre
    private Vector2[] initVertices;
    
    // The edges connecting the vertices
    private List<Pair<Vector2, Vector2>> edges;
        
    /**
     * Creates a new PhysicsBox for collision with the awesome car.
     * @param initMass     The mass of the box.
     * @param sideLength   The length of each side of the box.
     * @param initPosition The box's initial position.
     * @param initRotation The box's initial rotation.
     * @param initVel      The box's initial linear velocity.
     * @param initOmega    The box's initial rotational velocity.
     * @param elastic      Whether the box is elastic.
     */
    public PhysicsBox(float sideLength, Vector2 initPosition, float initRotation, float initMass, Vector2 initVel, float initOmega, bool elastic) : base(initMass, initPosition, initRotation, initVel, initOmega)
    {
        circleRadius = Mathf.Sqrt(((sideLength / 2) * (sideLength / 2)) +
                        ((sideLength / 2) * (sideLength / 2)));
        isElastic = elastic;
        sideLen = sideLength; 
        
        //  The box's vertices are arranged as below.
        //
        //  V0 ________V3
        //    |        |
        //    |        |
        //    |        |
        //    |________|
        //  V1         V2
        initVertices = new Vector2[4];
        initVertices[0] = new Vector2(-sideLength / 2, sideLength / 2);
        initVertices[1] = new Vector2(-sideLength / 2, -sideLength / 2);
        initVertices[2] = new Vector2(sideLength / 2, -sideLength / 2);
        initVertices[3] = new Vector2(sideLength / 2, sideLength / 2);
        
        vertices = new Vector2[4];
        edges = new List<Pair<Vector2, Vector2>>(4);
        
        for(int i = 0; i < initVertices.Length; i++)
        {
            vertices[i] = initVertices[i];
        }
        
        for(int i = 0; i < vertices.Length; i++)
        {
            Vector2 nextV = i == vertices.Length - 1 ? vertices[0] : vertices[i + 1];
            Pair<Vector2, Vector2> edge = new Pair<Vector2, Vector2>(vertices[i], nextV);
            edges.Add(edge);
        }
        updateVertices();
    }
    
    public override void update(float deltaTime)
    {
        Vector2 position = getPosition();
        rotation += angularVelocity * deltaTime * (float)(180/Mathf.PI);
        position.x += velocity.x * deltaTime;
        position.y += velocity.y * deltaTime;

        updateVertices();
    }

    public override void updateVertices()  
    {
        Vector2 position = getPosition();
        // Update all of the vertex positions based on the box's current position
        // and rotation
        for(int i = 0; i < vertices.Length; i++)
        {
            vertices[i].x = initVertices[i].x;
            vertices[i].y = initVertices[i].y;
            vertices[i].rotate(rotation);
            vertices[i].x += position.x;
            vertices[i].y += position.y;
        }
    }
    
    public override Vector2[] getVertices()
    {
        return vertices;
    }
    
    public override List<Pair<Vector2, Vector2>> getEdges()
    {
        return edges;
    }
    
    public override double getBoundingCircleRadius()
    {
        return circleRadius;
    }
    
    public override Vector2 getCentreOfMass()
    {
        return getPosition();
    }
    
    public override float getMomentOfInertia()
    {
        return (float)((this.mass * (Mathf.Pow(sideLen, 2) + Mathf.Pow(sideLen, 2))) / 12);
    }

    public override void setAngularVelocity(float vel)
    {
        
    }
}
