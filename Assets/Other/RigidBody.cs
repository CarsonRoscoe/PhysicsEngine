/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


using System.Collections.Generic;
using UnityEngine;
/**
*
* @author Shane
*/
public abstract class RigidBody {
    public float mass;
    public Vector2 position;
    public float rotation;
    public Vector2 velocity;
    public float angularVelocity;
    protected Vector2[] vertices;

    public RigidBody(float initMass, Vector2 initPosition, float initRotation,
                       Vector2 initVel, float initAngleVel)
    {
        mass = initMass;
        position = initPosition;
        rotation = initRotation;
        velocity = initVel;
        angularVelocity = initAngleVel;
    }
    
    /**
     * Gets the list of vertices for the given rigid object.
     * @return The list of vertices for this rigid object.
     */
    public abstract Vector2[] getVertices();
    
    /**
     * Gets the array of all edges.
     * @return The edges array.
     */
    public abstract List<Pair<Vector2, Vector2>> getEdges();
    
    /**
     * Gets the radius of the bounding circle to be used for broadphase collision
     * detection.
     * @return The radius of the object's bounding circle.
     */
    public abstract double getBoundingCircleRadius();
    
    /**
     * Update this object's physical state with the given time interval.
     * @param time The time interval over which to integrate.
     */
    public abstract void update(float time);
    
    /**
     * Update this object's vertex positions based on its current position. 
     */
    public abstract void updateVertices();
    
    /**
     * Gets the object's centre of mass.
     * @return The object's centre of mass.
     */
    public abstract Vector2 getCentreOfMass();
    
    /**
     * Gets this object's moment of inertia.
     * @return The object's moment of inertia.
     */
    public abstract float getMomentOfInertia();

    public abstract void setAngularVelocity(float vel);
    
    public Vector2 getPosition()
    {
        return position;
    }
    
    /**
     * Updates the position to the new position and updates the object's vertices.
     * @param newPos
     * @return 
     */
    public Vector2 setPosition(Vector2 newPos)
    {
        position = newPos;
        updateVertices();
        
        return position;
    }
}
