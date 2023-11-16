package empires.input;

import lmu.utils.EnumGetter;

public enum DSType {
    READS("reads"),
    DS("ds"),
    FRAGS("frags");

    String name;

    DSType(String name) {
        this.name = name;
        EnumGetter.add(name, this);
    }

    public static DSType get(String n) {
        return EnumGetter.get(n, true);
    }
}
