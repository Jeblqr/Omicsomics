"""Add pipeline_config to runs and custom_pipelines table

Revision ID: 0003_add_custom_pipelines
Revises: 0002_add_runs_and_datafiles
Create Date: 2025-11-09

"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSON


# revision identifiers, used by Alembic.
revision = "0003_add_custom_pipelines"
down_revision = "0002_add_runs_and_datafiles"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Add pipeline_config column to runs table
    op.add_column("runs", sa.Column("pipeline_config", JSON, nullable=True))

    # Create custom_pipelines table
    op.create_table(
        "custom_pipelines",
        sa.Column("id", sa.Integer(), autoincrement=True, nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("description", sa.Text(), nullable=False, server_default=""),
        sa.Column("definition", JSON, nullable=False),
        sa.Column(
            "category", sa.String(length=100), nullable=False, server_default="custom"
        ),
        sa.Column("is_public", sa.Boolean(), nullable=False, server_default="false"),
        sa.Column("owner_id", sa.Integer(), nullable=False),
        sa.Column("template_id", sa.String(length=100), nullable=True),
        sa.Column(
            "created_at",
            sa.DateTime(timezone=True),
            server_default=sa.text("now()"),
            nullable=False,
        ),
        sa.Column(
            "updated_at",
            sa.DateTime(timezone=True),
            server_default=sa.text("now()"),
            nullable=False,
        ),
        sa.ForeignKeyConstraint(
            ["owner_id"],
            ["users.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
    )
    op.create_index(
        op.f("ix_custom_pipelines_id"), "custom_pipelines", ["id"], unique=False
    )


def downgrade() -> None:
    op.drop_index(op.f("ix_custom_pipelines_id"), table_name="custom_pipelines")
    op.drop_table("custom_pipelines")
    op.drop_column("runs", "pipeline_config")
