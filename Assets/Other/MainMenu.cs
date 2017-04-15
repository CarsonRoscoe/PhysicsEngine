/*
public class MainMenu implements Screen {
    final Game game;

    OrthographicCamera camera;

    public MainMenu(final Game gam) {
        game = gam;

        camera = new OrthographicCamera();
        camera.setToOrtho(false, 800, 800);

    }

    @Override
    public void render(float deltaTime)
    {
        Gdx.gl.glClearColor(0, 0, 0.2f, 1);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);

        camera.update();
        game.batch.setProjectionMatrix(camera.combined);

        game.batch.begin();
        game.font.draw(game.batch, "Project Terminus", 100, 150);
        game.font.draw(game.batch, "Click to begin...", 100, 100);
        game.batch.end();

        if (Gdx.input.isTouched())
        {
            game.setScreen(new ProjectTerminus(game));
            dispose();
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