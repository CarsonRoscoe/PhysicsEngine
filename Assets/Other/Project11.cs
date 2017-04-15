
/*
public class Project11 //: Screen
{
    //final Game game;
    private OrthographicCamera camera;
    private ShapeRenderer shapeRenderer;
    private Vector2 initialVelocityU;
    private Vector2 initialVelocityV;
    private Vector2 initialCircle1Position;
    private Vector2 initialCircle2Position;
    private Vector2 normal;
    private Circle circle1;
    private Circle circle2;
    private float e; //coefficient of restitution
    private float impulse;
    private float pix;
    private float pfx;
    private float pfy;
    private float circle1KEi;
    private float circle1KEf;
    private float circle2KEi;
    private float circle2KEf;
    private float rect1RKEf;
    private float rect2RKEf;
    private float initialOrbital1;
    private float initialOrbital2;
    private float finalOrbital1;
    private float finalOrbital2;
    private float finalRotational1;
    private float finalRotational2;
    private int colisionCount;
    private Boolean running;
    private Boolean collided;
    private float rect1W;
    private float rect2W;
    private float rect1Rot;
    private float rect2Rot;

    public Project11(Game game)
    {
        this.game = game;

        camera = new OrthographicCamera();
        camera.setToOrtho(false, Gdx.graphics.getWidth(), Gdx.graphics.getHeight());
        camera.position.sub(Gdx.graphics.getWidth() / 2, Gdx.graphics.getHeight() / 2, 0);
        shapeRenderer = new ShapeRenderer();
        initialCircle1Position = new Vector2(-380, 0);
        initialCircle2Position = new Vector2(80, 0);
        circle1 = new Circle(initialCircle1Position, 1, Vector2.Zero, 40);
        circle2 = new Circle(initialCircle2Position, 1, Vector2.Zero, 40);

        setupInitialState();
    }

    private void setupInitialState()
    {
        collided = false;
        running = false;
        e = 0;
        impulse = 0;
        pix = 0;
        pfx = 0;
        pfy = 0;
        colisionCount = 0;
        initialVelocityU = new Vector2(100, 0);
        initialVelocityV = new Vector2(0, 0);
        circle1.position = new Vector2(initialCircle1Position);
        circle2.position = new Vector2(initialCircle2Position);
        circle1.velocity = new Vector2(initialVelocityU);
        circle2.velocity =  new Vector2(initialVelocityV);
        normal = new Vector2(0,0);
        circle1.mass = 1;
        circle2.mass = 1;

        circle1KEi = 0;
        circle1KEf = 0;
        circle2KEf = 0;
        circle2KEi = 0;

        rect1W = 0;
        rect2W = 0;
        rect1Rot = 0;
        rect2Rot = 0;
        rect1RKEf = 0;
        rect2RKEf = 0;
    }

    @Override
    public void render(float deltaTime)
    {
        Gdx.gl.glClearColor( 0.2f,  0.2f, 0.2f, 1);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);

        camera.update();
        game.batch.setProjectionMatrix(camera.combined);
        shapeRenderer.setProjectionMatrix(camera.combined);

        shapeRenderer.begin(ShapeRenderer.ShapeType.Line);
        shapeRenderer.setColor(1, 1, 0, 1);
        shapeRenderer.line(-Gdx.graphics.getWidth() / 2, 0, Gdx.graphics.getWidth() / 2, 0); // X Axis
        shapeRenderer.line(0, -Gdx.graphics.getHeight()/2, 0, Gdx.graphics.getHeight()/2); // Y Axis
        shapeRenderer.end();

        shapeRenderer.begin(ShapeRenderer.ShapeType.Filled);
        shapeRenderer.setColor(0, 1, 0, 1);
        shapeRenderer.translate(circle1.position.x, circle1.position.y, 0);
        shapeRenderer.rotate(0,0,1, rect1Rot);
        shapeRenderer.rect(-circle1.radius/2.0f, -circle1.radius/2.0f, circle1.radius, circle1.radius);
        shapeRenderer.rotate(0,0,1, -rect1Rot);
        shapeRenderer.translate(-circle1.position.x, -circle1.position.y, 0);
        shapeRenderer.setColor(0, 0, 1, 1);
        shapeRenderer.translate(circle2.position.x, circle2.position.y, 0);
        shapeRenderer.rotate(0, 0, 1, rect2Rot);
        shapeRenderer.rect(-circle2.radius/2.0f, -circle2.radius/2.0f, circle2.radius, circle2.radius);
        shapeRenderer.rotate(0,0,1, -rect2Rot);
        shapeRenderer.translate(-circle2.position.x, -circle2.position.y, 0);
        shapeRenderer.end();

        game.batch.begin();
        game.font.drawMultiLine(game.batch, "Mass 1: " + circle1.mass  + " Mass 2: " + circle2.mass + "\n" +
                "Initial: u = " + initialVelocityU + ", v = " + initialVelocityV + "\n" +
                "Final: u = " + circle1.velocity + ", v = " + circle2.velocity + "\n\n" +
                "e: " + e + "\n" +
                "Impulse: " + impulse + "\n" +
                "Normal: " + normal + "\n" +
                "Current Normal = " + "(" + (circle1.position.x - circle2.position.x) + ", " + (circle1.position.y - circle2.position.y) + ")" +
                " = " + "(" + circle1.position.x + ", " + circle1.position.y + ")" + " - " + "(" + circle2.position.x + ", " + circle2.position.y + ")" + "\n" +
                "Collision Count: " + colisionCount + "\n\n" +
                "p_ix = " + (circle1.mass * initialVelocityU.x) + " + " + (circle2.mass * initialVelocityV.x) + " = " + pix + "\n" +
                "p_fx = " + (circle1.mass * circle1.velocity.x) + " + " + (circle2.mass * circle2.velocity.x) + " = " +  pfx + "\n" +
                "p_fy = " + (circle1.mass * circle1.velocity.y) + " + " + (circle2.mass * circle2.velocity.y) + " = " +  pfy + "\n" +
                "L_i = " + initialOrbital1 + " + " + initialOrbital2 + " = " + (initialOrbital1 + initialOrbital2) + "\n" +
                "L_f = " + finalOrbital1 + " + " + finalRotational1 + "\n+ " + finalOrbital2 + "+ " + finalRotational2 + " = " + (finalOrbital1 + finalRotational1 + finalOrbital2 + finalRotational2) + "\n\n" +
                "KE_i = " + circle1KEi + " + " + circle2KEi + " = " + (circle1KEi + circle2KEi) + "\n" +
                "KE_f = " + circle1KEf + " + " + rect1RKEf + " + " + circle2KEf + " + " + rect2RKEf +  " = " + (circle1KEf + circle2KEf + rect1RKEf + rect2RKEf) + "\n", -Gdx.graphics.getWidth()/2, Gdx.graphics.getHeight()/2.2f);
        game.batch.end();

        if (running)
            updateSimulation(deltaTime);

        handleInput();
    }

    private void adjustCircleMass(Circle circle, Boolean increase)
    {
        float mass = circle.mass;

        if (increase) ++mass;
        else --mass;

        if (mass > 10) mass = 10;
        else if (mass < 0) mass = 0;

        circle.mass = mass;
    }

    private void adjustCircleVelocity(Circle circle, Boolean increase)
    {
        float velocity = circle.velocity.x;

        if (increase) velocity += 10;
        else velocity -= 10;

        if (circle.position.x == initialCircle1Position.x && circle.position.y == initialCircle1Position.y)
            initialVelocityU.x = velocity;
        else if (circle.position.x == initialCircle2Position.x && circle.position.y == initialCircle2Position.y)
            initialVelocityV.x = velocity;

        circle.velocity.x = velocity;
    }

    private void handleInput()
    {

        if (Gdx.input.isKeyJustPressed(Input.Keys.Q))
        {
            adjustCircleMass(circle1, true);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.W))
        {
            adjustCircleMass(circle1, false);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.E))
        {
            adjustCircleMass(circle2, true);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.R))
        {
            adjustCircleMass(circle2, false);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.A))
        {
            adjustCircleVelocity(circle1, true);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.S))
        {
            adjustCircleVelocity(circle1, false);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.D))
        {
            adjustCircleVelocity(circle2, true);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.F))
        {
            adjustCircleVelocity(circle2, false);
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.Z))
        {
            e += 0.1;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.X))
        {
            e -= 0.1f;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.SPACE))
        {
            running = !running;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.BACKSPACE))
        {
            setupInitialState();
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.DOWN))
        {
            circle1.position.y -= 10;
        }

        if (Gdx.input.isKeyJustPressed(Input.Keys.UP))
        {
            circle1.position.y += 10;
        }
    }

    private void updateSimulation(float deltaTime)
    {
        circle1KEi = (0.5f) * circle1.mass * (float)Math.pow(initialVelocityU.len(), 2);
        circle2KEi = (0.5f) * circle2.mass * (float)Math.pow(initialVelocityV.len(), 2);

        circle1.position.x += (circle1.velocity.x * deltaTime);
        circle2.position.x += (circle2.velocity.x * deltaTime);
        circle1.position.y += (circle1.velocity.y * deltaTime);
        circle2.position.y += (circle2.velocity.y * deltaTime);

        rect1Rot += (rect1W * deltaTime);
        rect2Rot += (rect2W * deltaTime);

        //float distance = Vector2.dst(circle1.position.x, circle1.position.y, circle2.position.x, circle2.position.y);
        // AABB
        Boolean collision = (circle1.position.x < circle2.position.x + circle2.radius &&
                circle1.position.x + circle1.radius > circle2.position.x &&
                circle1.position.y < circle2.position.y + circle2.radius &&
                circle1.radius + circle1.position.y > circle2.position.y);

        if (collision && collided)
        {
            colisionCount++;
        }

        if (collision && !collided) // if collision happened
        {
            circle1.position.x = circle2.position.x - 40;

            Vector3 unitNormal = new Vector3(1, 0, 0);

            Vector3 P = new Vector3(0,0,0);
            P.x = circle2.position.x - 20;
            P.y = circle1.position.y + ((circle2.position.y - circle1.position.y)/2.0f);

            Vector3 r1 = new Vector3((P.x - circle1.position.x), (P.y - circle1.position.y), 0);
            Vector3 r2 = new Vector3((P.x - circle2.position.x), (P.y - circle2.position.y), 0);

            initialOrbital1 = r1.y * circle1.mass * circle1.velocity.x;
            initialOrbital2 = r2.y * circle2.mass * circle2.velocity.x;

            float momentOfInertia1 = (circle1.mass * (float)(Math.pow(40, 2) + Math.pow(40, 2)))/12.0f;
            float momentOfInertia2 = (circle2.mass * (float)(Math.pow(40, 2) + Math.pow(40, 2)))/ 12.0f;

            float vR = circle1.velocity.x - circle2.velocity.x;

            // calculate impulse J
            float coefficient = -vR * (e + 1.0f);
            float inverseMasses = (1.0f/circle1.mass) + (1.0f/circle2.mass);

            Vector3 elementOne = new Vector3(r1);
            elementOne = elementOne.crs(unitNormal);
            elementOne.x /= momentOfInertia1;
            elementOne.y /= momentOfInertia1;
            elementOne.z /= momentOfInertia1;
            elementOne = elementOne.crs(r1);
            float componentOne = unitNormal.dot(elementOne); //Vector3.dot(unitNormal.x, unitNormal.y, unitNormal.z, elementOne.x, elementOne.y, elementOne.z);

            Vector3 elementTwo = new Vector3(r2);
            elementTwo = elementTwo.crs(unitNormal);
            elementTwo.x /= momentOfInertia2;
            elementTwo.y /= momentOfInertia2;
            elementTwo.z /= momentOfInertia2;
            elementTwo = elementTwo.crs(r2);
            float componentTwo = unitNormal.dot(elementTwo); //Vector3.dot(unitNormal.x, unitNormal.y, unitNormal.z, elementTwo.x, elementTwo.y, elementTwo.z);

            float denominator = inverseMasses + componentOne + componentTwo;

            float J = coefficient * (1.0f/denominator);

            float Uf = J/circle1.mass + circle1.velocity.x;
            float Vf = -J/circle2.mass + circle2.velocity.x;

            pix = (circle1.mass * initialVelocityU.x) + (circle2.mass * initialVelocityV.x);
            pfx = (circle1.mass * circle1.velocity.x) + (circle2.mass * circle2.velocity.x);
            pfy = (circle1.mass * circle1.velocity.y) + (circle2.mass * circle2.velocity.y);

            impulse = J;
            colisionCount++;

            circle1.velocity.x = Uf * unitNormal.x;
            circle1.velocity.y = Uf * unitNormal.y;
            circle2.velocity.x = Vf * unitNormal.x;
            circle2.velocity.y = Vf * unitNormal.y;

            Vector3 angularVelocity1 = new Vector3(unitNormal);
            angularVelocity1.x = angularVelocity1.x * J;
            angularVelocity1.y = angularVelocity1.y * J;
            angularVelocity1.z = angularVelocity1.z * J;
            angularVelocity1 = (new Vector3(r1)).crs(angularVelocity1);
            rect1W = angularVelocity1.z * (1.0f/momentOfInertia1) * (float)(180.0f/Math.PI);

            Vector3 angularVelocity2 = new Vector3(unitNormal);
            angularVelocity2.x = angularVelocity2.x * -J;
            angularVelocity2.y = angularVelocity2.y * -J;
            angularVelocity2.z = angularVelocity2.z * -J;
            angularVelocity2 = (new Vector3(r2)).crs(angularVelocity2);
            rect2W = angularVelocity2.z * (1.0f/momentOfInertia2) * (float)(180.0f/Math.PI);

            circle1KEf = (0.5f) * circle1.mass * (float)Math.pow(circle1.velocity.x, 2);
            circle2KEf = (0.5f) * circle2.mass * (float)Math.pow(circle2.velocity.x, 2);

            rect1RKEf = (0.5f) * momentOfInertia1 * (float)Math.pow(rect1W * (Math.PI/180.0f), 2);
            rect2RKEf = (0.5f) * momentOfInertia2 * (float)Math.pow(rect2W * (Math.PI/180.0f), 2);



            finalOrbital1 = r1.y * circle1.mass * circle1.velocity.len();
            finalOrbital2 = r2.y * circle2.mass * circle2.velocity.len();
            finalRotational1 = momentOfInertia1 * (float)(rect1W * (Math.PI/180.0f));
            finalRotational2 = momentOfInertia2 * (float)(rect2W * (Math.PI/180.0f));

            collided = true;
        }

    }

    private void calculateImpulse()
    {
        float vR = circle1.velocity.x - circle2.velocity.x;
        impulse = -vR * (e + 1) * ((circle1.mass * circle2.mass)/(circle1.mass + circle2.mass));
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