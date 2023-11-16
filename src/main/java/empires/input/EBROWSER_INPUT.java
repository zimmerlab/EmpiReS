package empires.input;

import java.io.File;

public enum EBROWSER_INPUT
{
    PDATA("p_data.txt"),
    EXPRS("exprs.txt"),
    FDATA("f_data.txt"),
    TRUES("trues.txt")
    ;

    String filename;

    EBROWSER_INPUT(String n)
    {
        this.filename = n;
    }

    public File get(File d)
    {
        return new File(d, filename);
    }
}
