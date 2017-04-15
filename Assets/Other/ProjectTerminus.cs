
/**
 * Created by Dan on 4/5/2015.
 */
 /*
public class ProjectTerminus implements Screen
{
    // For testing purposes
    private final Color eBoxUnmolestedColor = new Color(0.f, 0.f, 255.f, 1.f);
    private final Color eBoxHitColor = new Color(255.f, 0.f, 0.f, 1.f);
    private Color eBoxCurrentColor = eBoxUnmolestedColor;
    
    // Debug drawing shapes; get rid of these once it's fixed
    private ArrayList<Pair<Vector2, Vector2>> collidingEdges = new ArrayList<Pair<Vector2, Vector2>>();
    
    final Game game;
    private final OrthographicCamera camera;
    private final ShapeRenderer shapeRenderer;
    private final PhysicsRect car;
    private final PhysicsRect driver;
    private final PhysicsRect tank;
    private final PhysicsBox elasticBox;
    private final PhysicsBox inelasticBox;
    private final ArrayList<RigidBody> walls = new ArrayList<RigidBody>(4);
    private final ArrayList<RigidBody> bodies;
    private final Color COMcolour;
    private final Color arrowColour;
    private Vector2 arrowPosition;

    public ProjectTerminus(Game game)
    {
        this.game = game;

        shapeRenderer = new ShapeRenderer();

        camera = new OrthographicCamera();
        camera.setToOrtho(false, Gdx.graphics.getWidth(), Gdx.graphics.getHeight());
        camera.position.sub(Gdx.graphics.getWidth() / 2, Gdx.graphics.getHeight() / 2, 0);

        car = new PhysicsRect(-300, 0, 200, 100, new Color(0.5f,0.5f,0.5f,1f), 1000, 0);
        driver = new PhysicsRect(-240, 20, 40, 40, new Color(0.5f, 0f, 0f, 1f), 0, 0);
        tank = new PhysicsRect(-380, 0, 40, 80, new Color(0f, 0.5f, 0f, 1f), 0, 0);

        COMcolour = new Color(0,0,0.5f,1.0f);
        arrowColour = new Color(0.6f,0.6f,0f,1.0f);

        arrowPosition = new Vector2(0,0);

        car.addChild(driver);
        car.addChild(tank);
                
        elasticBox = new PhysicsBox(100.f, new Vector2(170, 0), 0, 50.f,
                                    new Vector2(), 0.f, true);

        inelasticBox = new PhysicsBox(100.f, new Vector2(0, 170), 0, 50.f,
                new Vector2(), 0.f, false);
        
        bodies = new ArrayList<RigidBody>();
        bodies.add(car);
        bodies.add(elasticBox);
        bodies.add(inelasticBox);
                
        // Create the walls of the stage and add them to the bodies ArrayList
        walls.add(new PhysicsBox(800, new Vector2(-800, 0), 0, 0, new Vector2(), 0, false));
        walls.add(new PhysicsBox(800, new Vector2(800, 0), 0, 0, new Vector2(), 0, false));
        walls.add(new PhysicsBox(800, new Vector2(0, -800), 0, 0, new Vector2(), 0, false));
        walls.add(new PhysicsBox(800, new Vector2(0, 800), 0, 0, new Vector2(), 0, false));
        
        bodies.addAll(walls);
    }

    
    //Detect and handle any collisions.
    
    private void handleCollisions()
    {
        ArrayList<Pair<RigidBody, RigidBody>> bPhaseResults = doBroadPhase(bodies);
        
        // Handle any actual reported collisions
        for(Pair<RigidBody, RigidBody> collPair : bPhaseResults)
        {            
            // Check if either object in the collision is a wall; if so, we only
            // want to check if they collided
            if(walls.contains(collPair.getRight()) || walls.contains(collPair.getLeft()))
            {
                RigidBody b1 = collPair.getLeft();
                RigidBody b2 = collPair.getRight();
                
                // Determine the object that's colliding with the wall and test for
                // collision
                RigidBody affected = walls.contains(b1) ? b2 : b1;
                RigidBody wall = affected == b1 ? b2 : b1;
                CollisionInfo ci = new CollisionInfo();
                if(doNarrowPhase(collPair, ci, false))
                {
                    // Determine if the body is moving in the same direction
                    // as the collision normal; if so, we need to reverse it to
                    // push in the right direction
                    float dot = affected.velocity.dot(ci.normal);
                    if(dot > 0)
                    {
                        ci.normal.scl(-1);
                    }
                    
                    // Push the affected body out of the wall in the direction of
                    // the collision normal
                    ci.normal.scl(ci.depth);
                    affected.position.add(ci.normal);
                    affected.updateVertices();
                    
                    if(affected instanceof PhysicsRect)
                    {
                        for(PhysicsRect child : ((PhysicsRect)affected).childRects)
                        {
                            child.position.add(ci.normal);
                        }
                    }
                    
                    // Could do the full check, but fuck that
                    if(Math.abs(wall.getPosition().x) > 10) affected.velocity.x *= -1;
                    else affected.velocity.y *= -1;
                }
                continue;
            }
            
            // Check if the two objects are actually colliding and respond appropriately if so
            CollisionInfo ci = new CollisionInfo();
            if(doNarrowPhase(collPair, ci, true))
            {
                float elasticity = 0;
                RigidBody b1 = collPair.getLeft();
                RigidBody b2 = collPair.getRight();

                // Use this incredible instanceof operator to determine the correct
                // elasticity
                if (b1 instanceof PhysicsBox && b2 instanceof PhysicsBox)
                {
                    elasticity += (((PhysicsBox) b1).isElastic) ? 1 : 0;
                    elasticity += (((PhysicsBox) b2).isElastic) ? 1 : 0;
                    elasticity /= 2;
                }
                else if (b1 instanceof PhysicsBox && b2 instanceof PhysicsRect)
                {
                    if (((PhysicsBox)b1).isElastic) elasticity = 1.0f;
                }
                else if (b1 instanceof PhysicsRect && b2 instanceof PhysicsBox)
                {
                    if (((PhysicsBox)b2).isElastic) elasticity = 1.0f;
                }

                // Apply the appropriate collision response
                collide(elasticity, ci);
                System.out.println("CollisionInfo: " + ci);
            }
        }
    }
    
    //
    //Does broadphase collision detection on the given list of of rigid bodies.
    //@param rBodies The list of bodies to check.
    //@return A list of pairs of rigid bodies that might be colliding.
    //
    private ArrayList<Pair<RigidBody, RigidBody>> doBroadPhase(ArrayList<RigidBody> rBodies)    
    {
        // Check every body against every other body
        // If their bounding circles overlap, add them to a list for closer inspection
        ArrayList<Pair<RigidBody, RigidBody>> possiblyColliding = new ArrayList<Pair<RigidBody, RigidBody>>();
        for(int b1Index = 0; b1Index < rBodies.size(); b1Index++)
        {
            RigidBody body1 = rBodies.get(b1Index);
            for(int b2Index = b1Index + 1; b2Index < rBodies.size(); b2Index++)
            {
                RigidBody body2 = rBodies.get(b2Index);
                double dist = body1.getPosition().dst(body2.getPosition());
                if(dist < (body1.getBoundingCircleRadius() + body2.getBoundingCircleRadius()))
                {
                    Pair<RigidBody, RigidBody> collPair = new Pair<RigidBody, RigidBody>(body1, body2);
                    possiblyColliding.add(collPair);
                }
            }
        }

        return possiblyColliding;
    }
    
    //
    //Checks whether the two given bodies are actually in collision.
    //
    //Note that this currently only actually checks rectangles; the algorithm
    //could be used to check other polygons by removing the axis optimisation.
    //
    //@param bodies             Two bodies that might be colliding.
    //@param ciOut              A CollisionInfo object to hold the collision information.
    //@param getCollisionPoints If this is false, only fills in the collision depth and normal.
    //
    //
    private boolean doNarrowPhase(Pair<RigidBody, RigidBody> bodies, CollisionInfo ciOut,
                                    boolean getCollisionPoints)
    {
        Vector2[] b1Vertices = bodies.getLeft().getVertices();
        Vector2[] b2Vertices = bodies.getRight().getVertices();

        // Use the separating axis theorem to check for collisions
        // We know that we're using boxes, so we can get away with checking
        // just 4 axes
        float x;
        Vector2 axes[] = new Vector2[4];
        // Calculate the edge vector and then find the normal (flip the slope and change the sign)
        axes[0] = new Vector2(b1Vertices[1]).sub(b1Vertices[0]);
        x = axes[0].x;
        axes[0].x = axes[0].y;
        axes[0].y = -x;
        axes[0].nor();

        axes[1] = new Vector2(b1Vertices[2]).sub(b1Vertices[1]);
        x = axes[1].x;
        axes[1].x = axes[1].y;
        axes[1].y = -x;
        axes[1].nor();

        axes[2] = new Vector2(b2Vertices[1]).sub(b2Vertices[0]);
        x = axes[2].x;
        axes[2].x = axes[2].y;
        axes[2].y = -x;
        axes[2].nor();

        axes[3] = new Vector2(b2Vertices[2]).sub(b2Vertices[1]);
        x = axes[3].x;
        axes[3].x = axes[3].y;
        axes[3].y = -x;
        axes[3].nor();

        int aIndex;
        int minAxisIdx = 0;
        float minTranslation = Float.MAX_VALUE;
        for(aIndex = 0; aIndex < axes.length; aIndex++)
        {
            Pair<Vector2, Vector2> minMax1 = new Pair<Vector2, Vector2>();
            Pair<Vector2, Vector2> minMax2 = new Pair<Vector2, Vector2>();

            Pair<Float, Float> proj1 = new Pair<Float, Float>();
            Pair<Float, Float> proj2 = new Pair<Float, Float>();

            getMinMax(bodies.getLeft(), axes[aIndex], minMax1, proj1);
            getMinMax(bodies.getRight(), axes[aIndex], minMax2, proj2);

            // If there's a gap between the projected vectors, then there was no collision
            // so return false
            if(proj1.getRight() < proj2.getLeft() || 
               proj2.getRight() < proj1.getLeft())
            {
                return false;
            }

            // There was no gap, so check how far the penetration is on this
            // axis
            float translation;
            translation = Math.min(proj1.getRight(), proj2.getRight()) -
                          Math.max(proj1.getLeft(), proj2.getLeft());                    
            if (translation < minTranslation)
            {
                minTranslation = translation;
                minAxisIdx = aIndex;
            }
        }
        
        // If we didn't specify to only detect collisions, fill in the given
        // CollisionInfo object
        if(getCollisionPoints)
        {
            // Flip the normal such that it always points from the first to the
            // second body in collision
            Vector2 obj1To2 = new Vector2(bodies.getRight().getPosition())
                                         .sub(bodies.getLeft().getPosition());
            if(new Vector2(obj1To2).nor().dot(axes[minAxisIdx]) < 0)
                    axes[minAxisIdx].scl(-1);

            // Determine the referent and incident faces
            Vector2 e1Max = new Vector2();
            Pair<Vector2, Vector2> bestEdge1 = getBestEdge(new Vector2(axes[minAxisIdx]),
                                                bodies.getLeft(), e1Max);
            Vector2 e2Max = new Vector2();
            Pair<Vector2, Vector2> bestEdge2 = getBestEdge(new Vector2(axes[minAxisIdx]).scl(-1),
                                                bodies.getRight(), e2Max);

            // Clip all of the points outside the overlapping range between the
            // two objects and return these as the collision points
            ArrayList<Vector2> collisionPoints = findCollisionPoints(axes[minAxisIdx],
                                                    bestEdge1, bestEdge2,
                                                    e1Max, e2Max);
            // Out objects are rectangles; if the collision manifold is empty and there are
            // no collision points in it, we've got an edge/edge collision, so add a point
            // halfway between objects 1 and 2's centres
            if(collisionPoints.isEmpty())
            {
                obj1To2.scl(0.5f);
                obj1To2.add(bodies.getLeft().position);
                collisionPoints.add(obj1To2);
            }
            
            ciOut.manifold = collisionPoints;
        }
        
        // Fill in the information determined in the check regardless of whether
        // we're getting the collision points.
        ciOut.depth = minTranslation;
        ciOut.normal = axes[minAxisIdx];
        ciOut.bodies = bodies;
        return true;
    }

    // clips the line segment points v1, v2
    // if they are past o along n
    ArrayList<Vector2> clip(Vector2 v1, Vector2 v2, Vector2 n, double o)
    {
        ArrayList<Vector2> cp = new ArrayList<Vector2>();
        double d1 = n.dot(v1) - o;
        double d2 = n.dot(v2) - o;
        // if either point is past o along n
        // then we can keep the point
        if (d1 >= 0.0) cp.add(v1);
        if (d2 >= 0.0) cp.add(v2);
        // finally we need to check if they
        // are on opposing sides so that we can
        // compute the correct point
        if (d1 * d2 < 0.0) {
            // if they are on different sides of the
            // offset, d1 and d2 will be a (+) * (-)
            // and will yield a (-) and therefore be
            // less than zero
            // get the vector for the edge we are clipping
            Vector2 e = new Vector2(v2.x - v1.x, v2.y - v1.y);
            // compute the location along e
            double u = d1 / (d1 - d2);
            e.x *= u;
            e.y *= u;
            e = e.add(v1);
            // add the point
            cp.add(e);
        }

        return cp;
    }

    private ArrayList<Vector2> findCollisionPoints(Vector2 normal,
                                                   Pair<Vector2, Vector2> edge1, Pair<Vector2, Vector2> edge2,
                                                   Vector2 e1MaxPoint, Vector2 e2MaxPoint)
    {
        Pair<Vector2, Vector2> refEdge, incEdge;
        Vector2 refMaxPoint;
        Vector2 edge1Vector = new Vector2(edge1.getRight().x - edge1.getLeft().x, edge1.getRight().y - edge1.getLeft().y);
        Vector2 edge2Vector = new Vector2(edge2.getRight().x - edge2.getLeft().x, edge2.getRight().y - edge2.getLeft().y);
        Boolean flip = false;

        // determine reference and incident edges
        if (Math.abs(edge1Vector.dot(normal)) <= Math.abs(edge2Vector.dot(normal)))
        {
            refEdge = edge1;
            incEdge = edge2;
            refMaxPoint = e1MaxPoint;
            System.out.println("no flip");
        }
        else
        {
            incEdge = edge1;
            refEdge = edge2;
            refMaxPoint = e2MaxPoint;
            flip = true;
            System.out.println("flip");
        }

        Vector2 refVector = new Vector2(refEdge.getRight().x - refEdge.getLeft().x, refEdge.getRight().y - refEdge.getLeft().y);
        refVector.nor();

        double offset1 = refVector.dot(refEdge.getLeft());
        // clip incident edge by the first vertex of the reference edge
        ArrayList<Vector2> clippedPoints = clip(incEdge.getLeft(), incEdge.getRight(), refVector, offset1);
        if (clippedPoints.size() < 2) return null; // failure

        double offset2 = refVector.dot(refEdge.getRight());
        // clip incident edge by the second vertex of reference edge in opposite direction
        clippedPoints = clip(clippedPoints.get(0), clippedPoints.get(1), new Vector2(-refVector.x, -refVector.y), -offset2);
        if (clippedPoints.size() < 2) return null; // failure

        // clip points past the reference edge along the reference edge's normal
        Vector2 refEdgeNormal = new Vector2(refVector.y, -refVector.x);
        if (flip) refEdgeNormal.scl(-1);

        double offset3 = refEdgeNormal.dot(refMaxPoint);
        ArrayList<Vector2> pointsToRemove = new ArrayList<Vector2>();

        for (Vector2 point : clippedPoints)
        {
            if (refEdgeNormal.dot(point) - offset3 < 0.0)
                pointsToRemove.add(point);
        }

        clippedPoints.removeAll(pointsToRemove);
        return clippedPoints;
    }
    
   //
   //Gets the edge most perpendicular to the given normal.
   //@param normal      The normal against which to compare the edge.
   //@param rect        The rectangle containing the edges.
   //@param furthestOut A Vector2 in which to store the best furthest point. 
   //@return The best edge.
   //
    private Pair<Vector2, Vector2> getBestEdge(Vector2 normal, RigidBody rect, Vector2 furthestOut)
    {
        Vector2[] vertices = rect.getVertices();
        
        // Find the vertex furthest along the normal
        float maxProj = vertices[0].dot(normal);
        int furthestIndex = 0;
        for(int i = 1; i < vertices.length; i++)
        {
            float dot = vertices[i].dot(normal);
            if(dot > maxProj)
            {
                maxProj = dot;
                furthestIndex = i;
            }
        }
        
        // Get the left and right edges containing this vertex
        Vector2 furthest = vertices[furthestIndex];
        Vector2 leftNeighbour = furthestIndex == 0 ? vertices[vertices.length - 1] : vertices[furthestIndex - 1];
        Vector2 rightNeighbour = furthestIndex == vertices.length - 1 ? vertices[0]: vertices[furthestIndex + 1];
        
        Vector2 leftEdge = (new Vector2(furthest).sub(leftNeighbour)).nor();
        Vector2 rightEdge = (new Vector2(furthest).sub(rightNeighbour)).nor();
        
        float lDot = leftEdge.dot(normal);
        float rDot = rightEdge.dot(normal);

        furthestOut.x = furthest.x;
        furthestOut.y = furthest.y;
        
        return lDot < rDot ? new Pair<Vector2, Vector2>(leftNeighbour, furthest)
                           : new Pair<Vector2, Vector2>(furthest, rightNeighbour);
    }
    
     //
     //Gets the minimum and maximum points projected on a given axis.
     //@param rect      The rectangle for which to get the points.
     //@param axis      The axis along which the rectangle's points will be projected.
     //@param pointsOut A pair in which to store the points. 
     //                 Left is the minimum point, right the maximum.
     //@param projectionsOut A pair in which to store the projections. Left is the
     //                      minimum projection, right the maximum.
     //
    private void getMinMax(RigidBody rect, Vector2 axis,
                                             Pair<Vector2, Vector2> pointsOut,
                                             Pair<Float, Float> projectionsOut)
    {
        final Vector2[] vertices = rect.getVertices();
        final int numVertices = vertices.length;
        Vector2 minPoint = vertices[0];
        Vector2 maxPoint = vertices[0];
        float minProj = vertices[0].dot(axis);
        float maxProj = vertices[0].dot(axis);
        
        // Project each point on to the axis of projection and determine if it
        // should be the new maximum or minimum point
        for(int i = 1; i < numVertices; i++)
        {
            float dot = vertices[i].dot(axis);
            if(dot < minProj)
            {
                minPoint = vertices[i];
                minProj = dot;
            }
            else if(dot > maxProj)
            {
                maxPoint = vertices[i];
                maxProj = dot;
            }
        }
        
        pointsOut.setLeft(minPoint);
        pointsOut.setRight(maxPoint);
        
        projectionsOut.setLeft(minProj);
        projectionsOut.setRight(maxProj);
    }
    
    //
    //Given an elasticity and collision information, appropriately model the
    //collision response.
    //@param elasticity    The elasticity of the collision.
    //@param collisionInfo An object containing the penetration depth, normal and
    //                     bodies involved in the collision.
    //
    private void collide(final float elasticity, final CollisionInfo collisionInfo)
    {
        // Check if we even have enough info to respond to the collision
        if (collisionInfo.manifold == null || collisionInfo.manifold.isEmpty())
        {
            return;
        }

        RigidBody b1 = collisionInfo.bodies.getLeft();
        RigidBody b2 = collisionInfo.bodies.getRight();
        
        // Separate the bodies
        // The algorithm always ensures that the normal points from body 1 to 2
        // We also add 5% to space them out a bit more
        Vector2 body1Push = new Vector2(collisionInfo.normal).scl(-0.505f);
        Vector2 body2Push = new Vector2(collisionInfo.normal).scl(0.505f);
        
        b1.position.add(body1Push);
        b2.position.add(body2Push);
        
        b1.updateVertices();
        b2.updateVertices();
        
        Vector3 unitNormal = new Vector3(-collisionInfo.normal.x, -collisionInfo.normal.y, 0);
        float b1Mass = b1.mass;
        float b2Mass = b2.mass;
        
        // If we're dealing with the car, handle its child rectangles appropriately
        if (b1 instanceof PhysicsRect)
        {
            b1Mass = ((PhysicsRect)b1).getTotalMass();
            for(PhysicsRect child : ((PhysicsRect)b1).childRects)
            {
                child.position.add(body1Push);
            }
        }
        else if (b2 instanceof PhysicsRect)
        {
            b2Mass = ((PhysicsRect)b2).getTotalMass();
            for(PhysicsRect child : ((PhysicsRect)b2).childRects)
            {
                child.position.add(body2Push);
            }
        }
        
        // Get the radial vectors
        Vector2 collPoint = collisionInfo.manifold.get(0);
        Vector3 p = new Vector3(collPoint.x, collPoint.y, 0);
        Vector3 r1 = new Vector3(b1.position.x, b1.position.y, 0).sub(p);
        Vector3 r2 = new Vector3(b2.position.x, b2.position.y, 0).sub(p);
        
        // Get some coefficients for the impulse calculation
        Vector3 normalVec3 = new Vector3(collisionInfo.normal.x, collisionInfo.normal.y, 0);
        float nDot1 = normalVec3.dot(new Vector3(r1).crs(normalVec3).scl(1 / b1.getMomentOfInertia()).crs(r1));
        float nDot2 = normalVec3.dot(new Vector3(r2).crs(normalVec3).scl(1 / b2.getMomentOfInertia()).crs(r2));
        float impulseDenominator = (1.f / b1Mass) + (1.f / b2Mass) + nDot1 + nDot2;

        // Get the initial velocities in the normal direction
        Vector2 v1N = new Vector2(collisionInfo.normal).scl(b1.velocity.dot(collisionInfo.normal));
        Vector2 v2N = new Vector2(collisionInfo.normal).scl(b2.velocity.dot(collisionInfo.normal));
       
        // Get the relative velocity and the impulse (as well as the impulse in the normal direction)
        Vector2 vRel = new Vector2(b1.velocity).sub(b2.velocity);
        Vector2 impulse = new Vector2(vRel).scl(-(elasticity + 1) / impulseDenominator);
        Vector2 impulseN = new Vector2(collisionInfo.normal).scl(impulse.dot(collisionInfo.normal));

        // Calculate the final linear velocities after the collision
        Vector2 vf1 = new Vector2(impulseN).scl(1 / b1Mass).add(v1N);
        Vector2 vf2 = new Vector2(impulseN).scl(-1 / b2Mass).add(v2N);
                
        b1.velocity = vf1;
        b2.velocity = vf2;
        
        // Calculate the final angular velocities
        float omegaDiff1 = (new Vector3(impulseN.x, impulseN.y, 0).crs(r1).z / b1.getMomentOfInertia());
        float omegaDiff2 = (new Vector3(impulseN.x, impulseN.y, 0).crs(r2).z / -b2.getMomentOfInertia());
        
        b1.angularVelocity = omegaDiff1;
        b2.angularVelocity = omegaDiff2;
    }
    
    @Override
    public void render(float deltaTime)
    {
        Gdx.gl.glClearColor( 0.2f,  0.2f, 0.2f, 1);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);

        for (RigidBody body : bodies)
        {
            body.update(deltaTime);
        }

        handleCollisions();
        arrowPosition = car.forcePosition;

        camera.update();
        game.batch.setProjectionMatrix(camera.combined);
        shapeRenderer.setProjectionMatrix(camera.combined);

        shapeRenderer.begin(ShapeRenderer.ShapeType.Line);
        shapeRenderer.setColor(1, 1, 0, 1);
        shapeRenderer.line(-Gdx.graphics.getWidth() / 2, 0, Gdx.graphics.getWidth() / 2, 0); // X Axis
        shapeRenderer.line(0, -Gdx.graphics.getHeight() / 2, 0, Gdx.graphics.getHeight() / 2); // Y Axis
        shapeRenderer.end();

        shapeRenderer.begin(ShapeRenderer.ShapeType.Filled); // draw car

        shapeRenderer.setColor(car.colour);
        shapeRenderer.identity();
        shapeRenderer.translate(car.getPosition().x, car.getPosition().y, 0.f);
        shapeRenderer.rotate(0, 0, 1, car.rotation);
        shapeRenderer.rect(-car.width / 2, -car.height / 2, car.width, car.height); // Car
        shapeRenderer.rotate(0, 0, 1, -car.rotation);
        shapeRenderer.translate(-car.getPosition().x, -car.getPosition().y, 0.f);

        shapeRenderer.setColor(driver.colour);
        shapeRenderer.translate(driver.getPosition().x, driver.getPosition().y, 0.f);
        shapeRenderer.rotate(0, 0, 1, driver.rotation);
        shapeRenderer.rect(-driver.width / 2, -driver.height / 2, driver.width, driver.height); // Driver
        shapeRenderer.rotate(0, 0, 1, -driver.rotation);
        shapeRenderer.translate(-driver.getPosition().x, -driver.getPosition().y, 0.f);

        shapeRenderer.setColor(tank.colour);
        shapeRenderer.translate(tank.getPosition().x, tank.getPosition().y, 0.f);
        shapeRenderer.rotate(0, 0, 1, tank.rotation);
        shapeRenderer.rect(-tank.width / 2, -tank.height / 2, tank.width, tank.height); // Driver
        shapeRenderer.rotate(0, 0, 1, -tank.rotation);
        shapeRenderer.translate(-tank.getPosition().x, -tank.getPosition().y, 0.f);

        shapeRenderer.setColor(COMcolour);
        shapeRenderer.circle(car.COM.x, car.COM.y, 4);
        shapeRenderer.end();
        // Draw the elastic box

        for (RigidBody body : bodies)
        {
            if (body instanceof PhysicsBox)
            {
                shapeRenderer.begin(ShapeRenderer.ShapeType.Filled);
                shapeRenderer.setColor(eBoxCurrentColor);
                PhysicsBox box = (PhysicsBox)body;
                shapeRenderer.translate(box.getPosition().x, box.getPosition().y, 0.f);
                shapeRenderer.rotate(0, 0, 1, box.rotation);
                shapeRenderer.rect(-box.sideLen / 2, -box.sideLen / 2,
                        box.sideLen, box.sideLen);
                shapeRenderer.rotate(0, 0, 1, -box.rotation);
                shapeRenderer.translate(-box.getPosition().x, -box.getPosition().y, 0.f);
                shapeRenderer.end();

                String elasticity = "1";

                if (!box.isElastic)
                {
                    elasticity = "0";
                }

                game.batch.begin();
                game.font.draw(game.batch, elasticity, box.getPosition().x, box.getPosition().y);
                game.batch.end();
            }
        }

        
        shapeRenderer.begin(ShapeRenderer.ShapeType.Filled); // draw car
        // Draw the force arrow if a force is currently being applied
        if (car.isForceOn(car.BACKWARD_FORCE))
        {
            arrowPosition = car.backForwardForceLocation();
            shapeRenderer.setColor(arrowColour);
            shapeRenderer.translate(arrowPosition.x, arrowPosition.y, 0.f);
            shapeRenderer.rotate(0, 0, 1, -(180 - car.rotation));
            shapeRenderer.rectLine(0, 0, 40, 0, 6);
            shapeRenderer.triangle(40, -10, 40, 10, 50, 0);
            shapeRenderer.rotate(0, 0, 1, (180 - car.rotation));
            shapeRenderer.translate(-arrowPosition.x, -arrowPosition.y, 0.f);
        }

        if (car.isForceOn(car.FORWARD_FORCE))
        {
            arrowPosition = car.backForwardForceLocation();
            shapeRenderer.setColor(arrowColour);
            shapeRenderer.translate(arrowPosition.x, arrowPosition.y, 0.f);
            shapeRenderer.rotate(0, 0, 1, car.rotation);
            shapeRenderer.rectLine(0, 0, 40, 0, 6);
            shapeRenderer.triangle(40, -10, 40, 10, 50, 0);
            shapeRenderer.rotate(0, 0, 1, -car.rotation);
            shapeRenderer.translate(-arrowPosition.x, -arrowPosition.y, 0.f);
        }

        if (car.isForceOn(car.TURNING_LEFT_FORCE))
        {
            arrowPosition = car.turnLeftForceLocation();
            shapeRenderer.setColor(arrowColour);
            shapeRenderer.translate(arrowPosition.x, arrowPosition.y, 0.f);
            shapeRenderer.rotate(0, 0, 1, car.rotation);
            shapeRenderer.rectLine(0, 0, 40, 0, 6);
            shapeRenderer.triangle(40, -10, 40, 10,  50, 0);
            shapeRenderer.rotate(0, 0, 1, -car.rotation);
            shapeRenderer.translate(-arrowPosition.x, -arrowPosition.y, 0.f);
        }

        if (car.isForceOn(car.TURNING_RIGHT_FORCE))
        {
            arrowPosition = car.turnRightForceLocation();
            shapeRenderer.setColor(arrowColour);
            shapeRenderer.translate(arrowPosition.x, arrowPosition.y, 0.f);
            shapeRenderer.rotate(0, 0, 1, car.rotation);
            shapeRenderer.rectLine(0, 0, 40, 0, 6);
            shapeRenderer.triangle(40, -10, 40, 10, 50, 0);
            shapeRenderer.rotate(0, 0, 1, -car.rotation);
            shapeRenderer.translate(-arrowPosition.x, -arrowPosition.y, 0.f);
        }

        shapeRenderer.end();
        
        for(Pair<Vector2, Vector2> edge: collidingEdges)
        {
            Vector2 p1 = edge.getLeft();
            Vector2 p2 = edge.getRight();
            shapeRenderer.begin(ShapeRenderer.ShapeType.Line);
            shapeRenderer.setColor(1, 1, 0, 1);
            shapeRenderer.line(p1.x, p1.y, p2.x, p2.y); // X Axis
            shapeRenderer.end();
        }
        
        game.batch.begin();

       //game.font.drawMultiLine(game.batch, "Car Position: " + car.getPosition() + "\n" +
       //        "Tank Position: " + tank.getPosition() + "\n" +
       //        "Driver Position: " + driver.getPosition() + "\n\n" +
       //        "COM Position: " + car.COM + "\n\n" +
       //        "Moment of Inertia: " + car.momentOfInertia + "\n" +
       //        "Car + Tank + Driver = " + "( " + car.mass + ", " + tank.mass + ", " + driver.mass + " )" + "\n\n" +
       //        "Force: " + car.force + "\n" +
       //        "v: " + car.velocity + "\n" +
       //        "a: " + car.linearAccel + "\n" +
       //        "theta: " +  (car.rotation * Math.PI/180) + "\n" +
       //        "omega: " + car.angularVelocity + "\n" +
       //        "alpha: " + car.angularAccel + "\n" +
       //        "Radial: " + car.radial + "\n", 100, Gdx.graphics.getHeight()/2.2f);

        game.font.draw(game.batch, "I'm a car.", car.getPosition().x, car.getPosition().y);
        game.batch.end();

        handleInput();
    }

    public void handleInput()
    {

        if (Gdx.input.isKeyJustPressed(Input.Keys.W))
        {
            car.mass += 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.S))
        {
            car.mass -= 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.Q))
        {
            tank.mass += 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.A))
        {
            tank.mass -= 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.E))
        {
            driver.mass += 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.D))
        {
            driver.mass -= 10;
        }

        if (Gdx.input.isKeyPressed(Input.Keys.UP))
        {
            car.centralForceOn(true);
        }
        else
        {
            car.centralForceOff(true);
        }

        if (Gdx.input.isKeyPressed(Input.Keys.DOWN))
        {
            car.centralForceOn(false);
        }
        else
        {
            car.centralForceOff(false);
        }

        if (Gdx.input.isKeyPressed(Input.Keys.RIGHT))
        {
            car.rightForceOn();
        }
        else
        {
            car.rightForceOff();
        }

        if (Gdx.input.isKeyPressed(Input.Keys.LEFT))
        {
            car.leftForceOn();
        }
        else
        {
            car.leftForceOff();
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.BACKSPACE))
        {
            car.reset(-300, 0, 200, 100, new Color(0.5f,0.5f,0.5f,1f), 1000, 0);
            driver.reset(-240, 20, 40, 40, new Color(0.5f, 0f, 0f, 1f), 0, 0);
            tank.reset(-380, 0, 40, 80, new Color(0f, 0.5f, 0f, 1f), 0, 0);
        }
    }

    @Override
    public void resize(int width, int height) {
    }

    @Override
    public void show() {
    }

    @Override
    public void hide() {
    }

    @Override
    public void pause() {
    }

    @Override
    public void resume() {
    }

    @Override
    public void dispose() {
    }
}
*/