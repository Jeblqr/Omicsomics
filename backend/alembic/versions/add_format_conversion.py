"""Add format conversion tables

Revision ID: add_format_conversion
Revises:
Create Date: 2025-01-10

"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "add_format_conversion"
down_revision = None  # Update this to point to the latest migration
branch_labels = None
depends_on = None


def upgrade():
    # Create format_conversions table
    op.create_table(
        "format_conversions",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("source_file_id", sa.Integer(), nullable=True),
        sa.Column("source_path", sa.String(), nullable=False),
        sa.Column("source_format", sa.String(), nullable=False),
        sa.Column("source_size_bytes", sa.Integer(), nullable=True),
        sa.Column("target_file_id", sa.Integer(), nullable=True),
        sa.Column("target_path", sa.String(), nullable=False),
        sa.Column("target_format", sa.String(), nullable=False),
        sa.Column("target_size_bytes", sa.Integer(), nullable=True),
        sa.Column(
            "conversion_path", postgresql.JSON(astext_type=sa.Text()), nullable=True
        ),
        sa.Column("conversion_mode", sa.String(), nullable=False),
        sa.Column("status", sa.String(), nullable=True, server_default="pending"),
        sa.Column("duration_seconds", sa.Float(), nullable=True),
        sa.Column("error_message", sa.String(), nullable=True),
        sa.Column("created_by", sa.Integer(), nullable=True),
        sa.Column(
            "created_at", sa.DateTime(), nullable=True, server_default=sa.text("now()")
        ),
        sa.Column("completed_at", sa.DateTime(), nullable=True),
        sa.Column("parameters", postgresql.JSON(astext_type=sa.Text()), nullable=True),
        sa.ForeignKeyConstraint(
            ["created_by"],
            ["users.id"],
        ),
        sa.ForeignKeyConstraint(
            ["source_file_id"],
            ["files.id"],
        ),
        sa.ForeignKeyConstraint(
            ["target_file_id"],
            ["files.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
    )
    op.create_index(
        op.f("ix_format_conversions_id"), "format_conversions", ["id"], unique=False
    )
    op.create_index(
        "ix_format_conversions_status", "format_conversions", ["status"], unique=False
    )
    op.create_index(
        "ix_format_conversions_created_by",
        "format_conversions",
        ["created_by"],
        unique=False,
    )

    # Create conversion_rules table
    op.create_table(
        "conversion_rules",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("from_format", sa.String(), nullable=False),
        sa.Column("to_format", sa.String(), nullable=False),
        sa.Column("method", sa.String(), nullable=False),
        sa.Column("script_path", sa.String(), nullable=True),
        sa.Column("command_template", sa.String(), nullable=True),
        sa.Column("avg_time_per_gb", sa.Float(), nullable=True),
        sa.Column("success_rate", sa.Float(), nullable=True, server_default="1.0"),
        sa.Column("usage_count", sa.Integer(), nullable=True, server_default="0"),
        sa.Column("requires_runtime", sa.String(), nullable=True),
        sa.Column(
            "requires_packages", postgresql.JSON(astext_type=sa.Text()), nullable=True
        ),
        sa.Column("is_active", sa.Boolean(), nullable=True, server_default="true"),
        sa.Column(
            "created_at", sa.DateTime(), nullable=True, server_default=sa.text("now()")
        ),
        sa.Column(
            "updated_at", sa.DateTime(), nullable=True, server_default=sa.text("now()")
        ),
        sa.Column("description", sa.String(), nullable=True),
        sa.Column("notes", sa.String(), nullable=True),
        sa.PrimaryKeyConstraint("id"),
    )
    op.create_index(
        op.f("ix_conversion_rules_id"), "conversion_rules", ["id"], unique=False
    )
    op.create_index(
        "ix_conversion_rules_from_to",
        "conversion_rules",
        ["from_format", "to_format"],
        unique=False,
    )


def downgrade():
    op.drop_index("ix_conversion_rules_from_to", table_name="conversion_rules")
    op.drop_index(op.f("ix_conversion_rules_id"), table_name="conversion_rules")
    op.drop_table("conversion_rules")

    op.drop_index("ix_format_conversions_created_by", table_name="format_conversions")
    op.drop_index("ix_format_conversions_status", table_name="format_conversions")
    op.drop_index(op.f("ix_format_conversions_id"), table_name="format_conversions")
    op.drop_table("format_conversions")
