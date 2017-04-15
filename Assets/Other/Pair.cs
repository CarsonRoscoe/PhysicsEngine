
/**
 *
 * @author Shane
 */
public class Pair<L, R> {
    private L left;
    private R right;
    
    public Pair(){}
    
    public Pair(L leftObj, R rightObj)
    {
        left = leftObj;
        right = rightObj;
    }
    
    public L getLeft() {return left;}
    public void setLeft(L newLeft) {left = newLeft;}
    
    public R getRight() {return right;}
    public void setRight(R newRight) {right = newRight;}
    
    public override string ToString()
    {
        return "Left: " + left.ToString() + ", Right: " + right.ToString();
    }
}
