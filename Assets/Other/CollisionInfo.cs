/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


using System.Collections.Generic;
using System.Text;
using UnityEngine;
/**
* Contains the information about a collision: the point and time of collision,
* objects after the objects have been moved back out of the collision, and the
* 
* @author Shane
*/
public class CollisionInfo {
    public List<Vector2> manifold;
    public float depth;
    public Vector2 normal;
    public Pair<RigidBody, RigidBody> bodies;
    
    /**
     * Default constructor.
     */
    public CollisionInfo(){}
    
    /**
     * Constructs a new CollisionInfo object with the specified information.
     * @param collidingBodies   The bodies in collision.
     * @param collisionManifold The collection of points in the collision.
     * @param PENETRATION_DEPTH The depth by which the bodies are overlapping.
     * @param normalDirection   The direction of the normal in the collision.
     */
    public CollisionInfo(Pair<RigidBody, RigidBody> collidingBodies, List<Vector2> collisionManifold, float PENETRATION_DEPTH, Vector2 normalDirection) {
        manifold = collisionManifold;
        normal = normalDirection;
        depth = PENETRATION_DEPTH;
        bodies = collidingBodies;
    }
    
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        sb.Append("Manifold: ");
        sb.Append(manifold);
        sb.Append("\nNormal: ");
        sb.Append(normal);
        sb.Append("\nDepth: ");
        sb.Append(depth);
        return sb.ToString();
    }
}
